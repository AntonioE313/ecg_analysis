import heartpy as hp
import scipy as sp
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.widgets import Button
import config
import debug_functions
import wfdb as wfdb
import neurokit2 as nk2
import math


class dataManager:
    def __init__(self):
        self.heartpy_params = {'wd': '',
                               'm': ''}
        print('DM - TEST')
        self.IO = self.IO()

        temp_dict = self.IO.wfdb_load()
        self.raw_data = temp_dict['data']

        self.data_wrapper = {'raw_data': self.raw_data,
                             'filtered_data': self.raw_data * 0}

        self.fs = temp_dict['fs']
        print('DM - self.fs= ', self.fs)
        beats_temp = self.preprocess_data()  # Placeholder for initializing beats wrapper

        # Dict for holding beat and data that is related index i to beat[i]
        # p_window_center_s - P window centre in samples
        self.beat_wrapper = {'beats': beats_temp,
                             'p_window_center_s': np.array([[0] * len(beats_temp)])[0],
                             'p_max_locs': np.array([[0] * len(beats_temp)])[0],
                             'premature': np.array([[0] * len(beats_temp)])[0],
                             'RR_intervals': np.array([[0.0] * len(beats_temp)])[0]}

        # Premature is 1 for premature 0 otherwise
        self.p_max_locs = np.array([[0] * len(self.beat_wrapper['beats'])])[
            0]  # Pre-allocate size beats[i] -> P_max_locs[i]

        self.SAECG_P_WINDOW_SIZE = int(config.p_window_duration * self.fs)  # Need to work out this number for even integer window size

        saecg_p_windows = np.tile(np.array(np.ones(shape=self.SAECG_P_WINDOW_SIZE)),
                                  (len(self.beat_wrapper['beats']), 1))

        saecg_p = np.array(np.ones(shape=self.SAECG_P_WINDOW_SIZE))

        # Wrapper for holding SAECG data
        self.saecg_wrapper = {'saecg_p': saecg_p,  # Signal averaged P wave window
                              'p_windows': saecg_p_windows}  # Windows used to calculate saecg_p

        # A dict holding the arrays corresponding to xcorr params, index i relates to beat[i]
        self.xcorr_params = {'le': np.array([[0] * len(self.beat_wrapper['beats'])])[0],
                             're': np.array([[0] * len(self.beat_wrapper['beats'])])[0],
                             'max_loc': np.array([[0] * len(self.beat_wrapper['beats'])])[0],
                             'max': np.array([[0.0] * len(self.beat_wrapper['beats'])])[0]}

        signal, waves = nk2.ecg_delineate(self.data_wrapper['raw_data'], self.heartpy_params['wd']['peaklist'],
                                          sampling_rate=3000, method="dwt", show=False,
                                          show_type='all')

        # Set all the Nan datapoints to 0 so they aren't a problem but don't mess up match with other beat by beat data
        for i in range(len(waves['ECG_P_Onsets'])):
            if math.isnan(waves['ECG_P_Onsets'][i]):
                waves['ECG_P_Onsets'][i] = 0
            if math.isnan(waves['ECG_P_Peaks'][i]):
                waves['ECG_P_Peaks'][i] = 0

        self.neurokit_params = {'signal': signal,
                                'waves': waves}

        temp_dict = self.find_p_peaks(self.beat_wrapper['beats'], self.beat_wrapper['p_window_center_s'])
        self.beat_wrapper['p_max_locs'] = temp_dict['p_max_locs']
        self.saecg_wrapper['p_windows'] = temp_dict['p_windows']
        self.compute_saecg()

    def preprocess_data(self):
        print('DM - preprocess_data')
        # run analysis
        self.data_wrapper['filtered_data'] = hp.filter_signal(self.data_wrapper['raw_data'], cutoff=0.05,
                                                              sample_rate=self.fs,
                                                              filtertype='notch')

        self.heartpy_params['wd'], self.heartpy_params['m'] = hp.process(hp.scale_data(self.data_wrapper['raw_data']),
                                                                         self.fs)

        # plt.figure()
        # plt.plot(self.data_wrapper['raw_data'])
        # plt.plot(self.data_wrapper['filtered_data'])
        # plt.show()

        print(self.heartpy_params['wd']['peaklist'])
        print(self.heartpy_params['wd']['peaklist'][1])

        print(type(self.heartpy_params['wd']['peaklist']))
        peaklist = self.heartpy_params['wd']['peaklist']
        beats = []

        for count, value in enumerate(peaklist, start=2):
            print(value, count)
            if count < len(peaklist):
                t1 = peaklist[count - 2]
                t2 = peaklist[count - 1]
                beats.append(self.data_wrapper['raw_data'][t1:t2])
                # beats.append(self.data_wrapper['filtered_data'][t1:t2])

        # display computed measures
        for measure in self.heartpy_params['m'].keys():
            print('%s: %f' % (measure, self.heartpy_params['m'][measure]))

        return beats

    # Need to make this function flexible not just for all beats, can take subset of beats
    def find_p_peaks(self, beats, p_window_center_s):
        # Prepare P windows array
        p_windows = np.tile(np.array(np.ones(shape=self.SAECG_P_WINDOW_SIZE-1)), (len(beats), 1))
        p_max_locs = np.array([[0] * len(beats)])[0]
        for i in range(len(beats)):
            # Generate the P windows
            beat_len = len(beats[i])

            if beat_len == 0:
                p_max_locs[i] = 0
                continue

            print('beat_len= ', beat_len)
            print('config.p_window_center= ', config.p_window_center)
            print('self.fs= ', self.fs)

            # p window center in samples (use as index)
            self.beat_wrapper['p_window_center_s'][i] = beat_len - int(config.p_window_center * self.fs)
            t1 = math.ceil(self.beat_wrapper['p_window_center_s'][i] - 0.5 * self.SAECG_P_WINDOW_SIZE)
            t2 = t1 + self.SAECG_P_WINDOW_SIZE
            debug_functions.array_inspect(p_windows, 'p_windows', False)
            # debug_functions.array_inspect(np.ndarray(beats), 'beats', False)
            print('t1 t2, ', t1, t2)
            print('p_window_center_s[i]', p_window_center_s[i])
            p_windows[i] = beats[i][t1:t2]
            p_max_locs_temp = sp.signal.find_peaks(p_windows[i], distance=len(p_windows[i]))[0]
            if len(p_max_locs_temp) > 1:  # If there are equal maxima with same value, take first one
                p_max_locs[i] = p_max_locs_temp[0]

            #   print('p_max_locs[i] = ', self.p_max_locs[i])
            p_max_locs[i] = p_max_locs[i] + self.beat_wrapper['p_window_center_s'][i] - round(
                0.5 * self.SAECG_P_WINDOW_SIZE)
            print('^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^')

        return {'self': self,
                'p_windows': p_windows,
                'p_max_locs': p_max_locs}

    def compute_saecg(self):
        for i in range(len(self.p_max_locs)):
            print()
            t1 = math.ceil(self.beat_wrapper['p_max_locs'][i] - 0.5 * self.SAECG_P_WINDOW_SIZE)
            t2 = t1 + self.SAECG_P_WINDOW_SIZE-1
            debug_functions.array_inspect(self.saecg_wrapper['p_windows'][i], 'self.saecg_wrapper[p_windows][i] ', False)
            debug_functions.array_inspect(self.beat_wrapper['beats'][i][t1:t2], 'self.beat_wrapper[beats][i][t1:t2]', False)
            self.saecg_wrapper['p_windows'][i] = self.beat_wrapper['beats'][i][t1:t2]
        self.saecg_wrapper['saecg_p'] = np.sum(self.saecg_wrapper['p_windows'], axis=0)
        self.saecg_wrapper['saecg_p'] = np.true_divide(self.saecg_wrapper['saecg_p'], len(self.beat_wrapper['beats']))
        debug_functions.array_inspect(self.saecg_wrapper['saecg_p'], 'saecg_p', False)
        return self

    # Compute the xcorr of SAECG with the array of beats using le and re of P template (specified by user)
    def compute_saecg_xcorr(self, le, re):
        saecg_local = self.saecg_wrapper['saecg_p'][int(le * self.fs):int(re * self.fs)]  # saecg cropped with le and re
        for i in range(len(self.beat_wrapper['beats'])):
            # percent change data
            pc_change_saecg_local = 100 * np.true_divide(np.diff(saecg_local), saecg_local[1:])
            pc_change_beats_i = 100 * np.true_divide(np.diff(self.beat_wrapper['beats'][i]),
                                                     self.beat_wrapper['beats'][i][1:])

            # Define P window edges as ratio of full beat length
            xcorr_p_window_le = 0.5
            xcorr_p_window_re = 0.9

            xcorr = np.correlate(pc_change_saecg_local, pc_change_beats_i, mode='same')
            xcorr_p_window = xcorr[round(xcorr_p_window_le * len(xcorr)):round(xcorr_p_window_re * len(xcorr))]

            # Set xcorr_params
            self.xcorr_params['max_loc'][i] = xcorr_p_window.tolist().index(np.amax(xcorr_p_window))

            self.xcorr_params['le'][i] = sp.signal.find_peaks(-xcorr_p_window[0:self.xcorr_params['max_loc'][i]],
                                                              distance=self.fs)[0]

            self.xcorr_params['re'][i] = sp.signal.find_peaks(-xcorr_p_window[self.xcorr_params['max_loc'][i]:],
                                                              distance=self.fs)[0] + self.xcorr_params['max_loc'][i]

            self.xcorr_params['max'][i] = xcorr[self.xcorr_params['max_loc'][i]]

            print('DM - self.xcorr_params[max][i]=', self.xcorr_params['max'][i])
            print('DM - self.xcorr_params[le][i]', self.xcorr_params['le'][i])
            print('DM - self.xcorr_params[re][i]', self.xcorr_params['re'][i])

            # Debug plots
            # plt.figure()
            # plt.plot(xcorr_p_window)
            # plt.plot(xcorr_p_window[0:self.xcorr_params['max_loc'][i]])
            # plt.plot(self.xcorr_params['le'][i], xcorr_p_window[self.xcorr_params['le'][i]], 'x')
            # plt.plot(self.xcorr_params['max_loc'][i], xcorr_p_window[self.xcorr_params['max_loc'][i]], 'o')
            # plt.show()

            p_window_le = round(0.5 * len(xcorr))  # loc of p window le in samples from start of beat
            print('DM - p_window_le', p_window_le)
            print('DM - 0.5*len(saecg_local)', 0.5 * len(saecg_local))
            print('DM - self.p_window_center_s[i]', self.beat_wrapper['p_window_center_s'][i])
            print('DM - p_window_max_loc=', self.xcorr_params['max_loc'][i])
            print('len(xcorr)', len(xcorr))
            print('len(self.beat_wrapper[beats][i])', len(self.beat_wrapper['beats'][i]))

            # Set the loc params to be in terms of beat rather than p window
            self.xcorr_params['max_loc'][i] = self.xcorr_params['max_loc'][i] + p_window_le
            self.xcorr_params['le'][i] = self.xcorr_params['le'][i] + p_window_le
            self.xcorr_params['re'][i] = self.xcorr_params['re'][i] + p_window_le

            # print('temp.tolist().index(np.amax(temp)) = ', temp.tolist().index(np.amax(temp)))

        print('DM - self.self.xcorr_params[max_loc][i] : ', self.xcorr_params['max_loc'][i])
        print('DM - self.xcorr_params[max][i] : ', self.xcorr_params['max_loc'][i])
        return self

    def premature_analysis(self, event):
        # Calculate RR intervals
        for i in range(len(self.beat_wrapper['beats'])):
            self.beat_wrapper['RR_intervals'][i] = len(self.beat_wrapper['beats'][i])

        # If the RR interval for beat i is < premature_signal_coefficient* the RR interval for beat (i-1) then note beat i as premature
        print('DM - Premature Analysis - self.beat_wrapper[RR_intervals]= ', self.beat_wrapper['RR_intervals'])
        for i in range(1, len(self.beat_wrapper['RR_intervals'])):
            if self.beat_wrapper['RR_intervals'][i] <= \
                    config.premature_signal_coefficient * self.beat_wrapper['RR_intervals'][i - 1]:
                self.beat_wrapper['premature'][i] = 1
            else:
                self.beat_wrapper['premature'][i] = 0

        NUMBER_OF_STANDARD_BEATS = len(self.beat_wrapper['beats']) - np.sum(self.beat_wrapper['premature'])

        # print('DM - premature beats= ', np.sum(self.beat_wrapper['premature']))
        # print('DM - standard beats= ', NUMBER_OF_STANDARD_BEATS)

        # Compute SAECG of non premature (standard) beats
        standard_beats = []
        for i in range(len(self.beat_wrapper['beats'])):
            if self.beat_wrapper['premature'][i] == 0:
                standard_beats.append(self.beat_wrapper['beats'][i])

        temp_dict = self.find_p_peaks(standard_beats, self.beat_wrapper['p_window_center_s'])
        standard_p_max_locs = temp_dict['p_max_locs']
        # (STATE VARIABLES FOR P PEAK LOCATION)
        # print('DM - standard_p_max_locs= ', standard_p_max_locs)
        standard_saecg_p_windows = []
        # print('DM - beat len 0 1 = ', len(self.beat_wrapper['beats'][0]), len(self.beat_wrapper['beats'][1]))
        for i in range(NUMBER_OF_STANDARD_BEATS):
            #    print('DM - i=', i)
            t1 = int(standard_p_max_locs[i] - 0.5 * self.SAECG_P_WINDOW_SIZE)
            t2 = int(standard_p_max_locs[i] + 0.5 * self.SAECG_P_WINDOW_SIZE)
            #    print('DM - len(standard_beats[i])= ', len(standard_beats[i]))
            #    print('DM - t1 t2 ', t1, t2)
            #    print('standard_beats[i][t1:t2]', standard_beats[i][t1:t2])
            #    print('standard_saecg_p_windows= ', standard_saecg_p_windows)
            standard_saecg_p_windows.append(standard_beats[i][t1:t2])
        standard_saecg_p_windows = np.sum(standard_saecg_p_windows, axis=0)
        standard_saecg_p_windows = np.true_divide(standard_saecg_p_windows, NUMBER_OF_STANDARD_BEATS)

        # Get percent change RR interval compared to beat before
        RRpc = np.array([[0.0] * len(self.beat_wrapper['RR_intervals'])])[0]
        RRpc[0] = 1
        for i in range(1, len(RRpc)):
            RRpc[i] = 100 * (self.beat_wrapper['RR_intervals'][i - 1] - self.beat_wrapper['RR_intervals'][i]) / \
                      self.beat_wrapper['RR_intervals'][i - 1]
        plt.figure()
        plt.title('Instantaneous % change in RR interval vs beat number')
        plt.plot(RRpc)
        plt.show()
        # Get interquartile range
        # iqr =

        return standard_saecg_p_windows

    def get_saecg_p(self):
        return self.saecg_p

    def get_p_max_locs(self):
        return self.p_max_locs

    def set_p_max_locs(self, p_max_locs):
        self.p_max_locs = p_max_locs

    class IO:
        def __init__(self):
            print('IO - init')

        # Save the passed dict as a csv file
        def dict_as_csv(self, d):
            pd.DataFrame(d).to_csv('test.csv')

        def csv_load(self, fs):
            print('IO - load data')
            signal_duration = config.signal_duration

            print('IO - config.data_filepath[-3:]=', config.data_filepath[-3:])
            if config.data_filepath[-3:] == 'csv':
                print('IO - csv load executed')
                raw_data = np.loadtxt(open(config.data_filepath, "rb"), skiprows=1)
                return raw_data[0:signal_duration * fs, 1]

        def wfdb_load(self):
            print('IO - wfdb load executed')
            record = wfdb.rdrecord(config.data_filepath)
            raw_data = record.p_signal[0:record.fs * config.signal_duration, 1]
            print('IO - record.fs= ', record.fs)

            # Debug plots
            # plt.figure()
            # plt.plot(raw_data)
            # plt.show()

            print('IO - type(raw data)=', type(raw_data))
            print('IO - np.shape(raw data)=', np.shape(raw_data))
            resampled_data = sp.signal.resample(raw_data, config.upsample_coefficient * len(raw_data))

            return {'data': resampled_data,
                    'fs': config.upsample_coefficient * record.fs}


class UIManager:
    def __init__(self, index, dm):
        self.current_index = index
        self.dm = dm
        self.raw_data = dm.data_wrapper['raw_data']
        self.p_max_locs = dm.p_max_locs
        self.fs = dm.fs
        self.mode = 1  # view mode for beat viewer

        self.fig, self.ax = plt.subplots(nrows=5, ncols=1, gridspec_kw={'height_ratios': [1, 1, 0.1, 0.1, 1]})
        self.fig.tight_layout()  # Configure the layout of UI elements

        self.axprev = plt.axes([0.7, 0.01, 0.08, 0.03])
        self.axnext = plt.axes([0.81, 0.01, 0.08, 0.03])
        self.axanalyze = plt.axes([0.05, 0.01, 0.08, 0.03])
        self.axmode = plt.axes([0.15, 0.01, 0.08, 0.03])
        self.axpremature = plt.axes([0.25, 0.01, 0.1, 0.03])

        self.bnext = Button(self.axnext, 'Next')
        self.bnext.on_clicked(lambda x: self.next_button_pushed(x))

        self.bprev = Button(self.axprev, 'Previous')
        self.bprev.on_clicked(lambda x: self.prev_button_pushed(x))

        self.banalyze = Button(self.axanalyze, 'Analyze')
        self.banalyze.on_clicked(lambda x: self.analyze_button_pushed(x))

        self.bmode = Button(self.axmode, 'Mode')
        self.bmode.on_clicked(lambda x: self.mode_button_pushed(x))

        self.bpremature = Button(self.axpremature, 'Premature Beat Analysis')
        self.bpremature.on_clicked(lambda x: self.premature_button_pushed(x))

        plt.axes(self.ax[0])
        self.ax[0].set_title('Raw Data')
        plt.plot(np.true_divide(range(len(sp.signal.resample(self.raw_data, len(self.raw_data)))), self.fs),
                 self.raw_data)

        plt.scatter(np.true_divide(self.dm.heartpy_params['wd']['peaklist'], self.fs),
                    self.raw_data[self.dm.heartpy_params['wd']['peaklist']], c='g')

        plt.scatter(np.true_divide(self.dm.neurokit_params['waves']['ECG_P_Peaks'], self.fs),
                    self.raw_data[self.dm.neurokit_params['waves']['ECG_P_Peaks']], marker='x', c='r')

        plt.axes(self.ax[1])
        title_string = 'Signal Averaged P Wave'
        self.ax[1].set_title(title_string)
        plt.plot(np.true_divide(range(len(self.dm.saecg_wrapper['saecg_p'])), self.fs),
                 self.dm.saecg_wrapper['saecg_p'])

        # Make the slider for left edge
        plt.axes(self.ax[2])
        self.slider_start = plt.Slider(self.ax[2], 'Start', 0,
                                       len(self.dm.saecg_wrapper['saecg_p']) / self.fs,
                                       valinit=0.1 * len(self.dm.saecg_wrapper['saecg_p']) / self.fs,
                                       color='blue')
        self.slider_start.on_changed(lambda x: self.slider_updated(x))
        # Make the slider for right edge
        plt.axes(self.ax[3])
        self.slider_end = plt.Slider(self.ax[3], 'End', 0, len(self.dm.saecg_wrapper['saecg_p']) / self.fs,
                                     valinit=0.9 * len(self.dm.saecg_wrapper['saecg_p']) / self.fs,
                                     color='blue')
        self.slider_end.on_changed(lambda x: self.slider_updated(x))

        plt.axes(self.ax[1])
        self.ax[1].vlines(self.slider_start.val,
                          self.dm.saecg_wrapper['saecg_p'][int(self.slider_start.val * self.fs)] - 0.02,
                          self.dm.saecg_wrapper['saecg_p'][int(self.slider_start.val * self.fs)] + 0.02, color='red')
        self.ax[1].vlines(self.slider_end.val,
                          self.dm.saecg_wrapper['saecg_p'][int(self.slider_end.val * self.fs)] - 0.02,
                          self.dm.saecg_wrapper['saecg_p'][int(self.slider_end.val * self.fs)] + 0.02, color='red')

        plt.axes(self.ax[4])
        self.ax[4].clear()
        title_string = 'Beat %d out of %d beats' % (self.get_current_index() + 1, len(dm.beat_wrapper['beats']))
        self.ax[4].set_title(title_string)

        plt.plot(np.true_divide(range(len(self.dm.beat_wrapper['beats'][self.current_index])),
                                self.fs),
                 self.dm.beat_wrapper['beats'][self.current_index])

        plt.plot(self.p_max_locs[self.get_current_index()] / self.fs,
                 self.dm.beat_wrapper['beats'][self.get_current_index()][
                     round(self.p_max_locs[self.get_current_index()])],
                 marker='x')

        plt.plot(self.dm.xcorr_params['max_loc'][self.get_current_index()] / self.fs,
                 self.dm.beat_wrapper['beats'][self.current_index][
                     self.dm.xcorr_params['max_loc'][self.get_current_index()]],
                 'o')

        mng = plt.get_current_fig_manager()
        mng.window.state("zoomed")
        plt.show()

    def refresh_saecg_view(self, saecg, title_string):
        print('UI - refresh_saecg_view')
        plt.axes(self.ax[1])
        self.ax[1].set_title(title_string)
        plt.plot(np.true_divide(range(len(saecg)), self.fs), saecg)
        plt.show()

    def refresh_beat_viewer(self):
        if self.get_current_index() >= 0:
            plt.axes(self.ax[4])
            self.ax[4].clear()
            title_string = 'Beat %d out of %d beats' % (self.get_current_index() + 1, len(dm.beat_wrapper['beats']))
            self.ax[4].set_title(title_string)

            if self.mode == 1:
                plt.plot(np.true_divide(range(len(self.dm.beat_wrapper['beats'][self.current_index])),
                                        self.fs),
                         self.dm.beat_wrapper['beats'][self.current_index])

                plt.plot(self.p_max_locs[self.get_current_index()] / self.fs,
                         self.dm.beat_wrapper['beats'][self.get_current_index()][
                             round(self.p_max_locs[self.get_current_index()])],
                         marker='x')

                plt.plot(
                    self.dm.xcorr_params['max_loc'][self.get_current_index()] / self.fs,
                    self.dm.beat_wrapper['beats'][self.current_index][
                        self.dm.xcorr_params['max_loc'][self.get_current_index()]],
                    'o')

            elif self.mode == 2:
                # Create the array which is the 3 beats index-1 through index + 1
                display_data = np.append(np.append(self.dm.beat_wrapper['beats'][self.current_index - 1],
                                                   self.dm.beat_wrapper['beats'][self.current_index]),
                                         self.dm.beat_wrapper['beats'][self.current_index + 1])
                plt.plot(np.true_divide(range(len(display_data)), self.fs), display_data)

                plt.plot((len(self.dm.beat_wrapper['beats'][self.current_index - 1]) + self.dm.xcorr_params['max_loc'][
                    self.get_current_index()]) / self.fs,
                         display_data[(len(self.dm.beat_wrapper['beats'][self.current_index - 1]) +
                                       self.dm.xcorr_params['max_loc'][self.get_current_index()])],
                         marker='o')

                plt.plot(
                    (len(self.dm.beat_wrapper['beats'][self.current_index - 1]) + self.dm.beat_wrapper['p_max_locs'][
                        self.current_index]) / self.fs,
                    display_data[
                        len(self.dm.beat_wrapper['beats'][self.current_index - 1]) + self.dm.beat_wrapper['p_max_locs'][
                            self.current_index]],
                    marker='x')

                plt.plot((len(self.dm.beat_wrapper['beats'][self.current_index - 1]) + self.dm.xcorr_params['re'][
                    self.current_index]) / self.fs,
                         self.dm.beat_wrapper['beats'][self.current_index][
                             self.dm.xcorr_params['re'][self.current_index]],
                         marker='*')

            plt.show()

    def next_button_pushed(self, event):
        # print('data type of index : ', type(self.get_current_index()))
        print('self.dm.xcorr_params[max_loc] = ', self.dm.xcorr_params['max_loc'])
        if self.get_current_index() < len(self.dm.beat_wrapper['beats']) - 1:
            plt.axes(self.ax[4])
            self.set_current_index(self.get_current_index() + 1)
            self.ax[4].clear()
            title_string = 'Beat %d out of %d beats' % (
                self.get_current_index() + 1, len(self.dm.beat_wrapper['beats']))
            self.ax[4].set_title(title_string)
            self.refresh_beat_viewer()

            print('Type of self.p_max_locs[self.get_current_index()]: ',
                  type(self.p_max_locs[self.get_current_index()]))
            print('Type of self.get_current_index() : ', type(self.get_current_index()))
            print('data type of index : ', type(self.get_current_index()))
            plt.show()

    def prev_button_pushed(self, event):
        if self.get_current_index() > 0:
            plt.axes(self.ax[4])
            self.set_current_index(self.current_index - 1)
            self.ax[4].clear()
            title_string = 'Beat %d out of %d beats' % (self.current_index + 1, len(dm.beat_wrapper['beats']))
            self.ax[4].set_title(title_string)
            self.refresh_beat_viewer()

    def analyze_button_pushed(self, event):
        self.dm.find_p_peaks(beats=self.dm.beat_wrapper['beats'],
                             p_window_center_s=self.dm.beat_wrapper['p_window_center_s'])
        self.dm.compute_saecg()
        self.dm.compute_saecg_xcorr(self.slider_start.val, self.slider_end.val)
        self.dm.IO.dict_as_csv(self.dm.xcorr_params)
        print('type(self)', type(self))
        print('UI - self.dm.xcorr_params[max_loc]=', self.dm.xcorr_params['max_loc'])
        print('UI - self.dm.xcorr_params[max][self.current_index]=', self.dm.xcorr_params['max'][self.current_index])

    def mode_button_pushed(self, event):
        if self.mode == 1:
            self.mode = 2
        elif self.mode == 2:
            self.mode = 1
        print('self.mode', self.mode)
        self.refresh_beat_viewer()

    def premature_button_pushed(self, x):
        print('UI - premature')
        standard_beats = self.dm.premature_analysis(x)
        self.refresh_saecg_view(standard_beats, 'SAECG of standard P-waves')

    def slider_updated(self, event):
        plt.axes(self.ax[1])
        self.ax[1].clear()
        title_string = 'Signal Averaged P Wave'
        self.ax[1].set_title(title_string)
        plt.plot(np.true_divide(range(len(self.dm.saecg_wrapper['saecg_p'])), self.fs),
                 self.dm.saecg_wrapper['saecg_p'])
        self.ax[1].vlines(self.slider_start.val,
                          self.dm.saecg_wrapper['saecg_p'][int(self.slider_start.val * self.fs)] - 0.02,
                          self.dm.saecg_wrapper['saecg_p'][int(self.slider_start.val * self.fs)] + 0.02, color='red')
        self.ax[1].vlines(self.slider_end.val,
                          self.dm.saecg_wrapper['saecg_p'][int(self.slider_end.val * self.fs)] - 0.02,
                          self.dm.saecg_wrapper['saecg_p'][int(self.slider_end.val * self.fs)] + 0.02, color='red')

    def get_current_index(self):
        return self.current_index

    def set_current_index(self, index):
        self.current_index = index


# START HERE
dm = dataManager()
ui = UIManager(1, dm)

plt.show()
