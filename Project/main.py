import heartpy as hp
import scipy as sp
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.widgets import Button



class dataManager:
    # ryan test 
    def __init__(self):
        self.heartpy_params = {'wd': '',
                               'm': ''}
        print('DM - TEST')

        self.IO = self.IO()

        self.fs = 250
        self.raw_data = self.IO.load_data(self.fs)
        self.beats = []
        self.preprocess_data()

        self.saecg_window_size = int(0.2 * self.fs)  # Watch for rounding error here, need to convert to int
        self.p_window_center = 0.15  # P window centre as in seconds back from end of beat
        self.p_window_center_s = np.array([[0] * len(self.beats)])[
            0]  # loc of p window centre in samples from start of beat
        self.p_max_locs = np.array([[0] * len(self.beats)])[0]  # Pre-allocate size beats[i] -> P_max_locs[i]
        self.saecg_p_windows = np.tile(np.array(np.ones(shape=(self.saecg_window_size))), (len(self.beats), 1))
        self.saecg_p = np.array(np.ones(shape=self.saecg_window_size))

        # A dict holding the arrays corresponding to xcorr params, index i relates to beat[i]
        self.xcorr_params = {'le': np.array([[0] * len(self.beats)])[0],
                             're': np.array([[0] * len(self.beats)])[0],
                             'max_loc': np.array([[0] * len(self.beats)])[0],
                             'max': np.array([[0.0] * len(self.beats)])[0]}

        self.find_p_peaks()
        self.compute_saecg()

    def preprocess_data(self):
        print('DM - preprocess_data')
        # run analysis
        self.heartpy_params['wd'], self.heartpy_params['m'] = hp.process(self.raw_data, self.fs)
        print(self.heartpy_params['wd']['peaklist'])
        print(self.heartpy_params['wd']['peaklist'][1])

        print(type(self.heartpy_params['wd']['peaklist']))
        peaklist = self.heartpy_params['wd']['peaklist']
        for count, value in enumerate(peaklist, start=2):
            print(value, count)
            if count < len(peaklist):
                t1 = peaklist[count - 2]
                t2 = peaklist[count - 1]
                self.beats.append(self.raw_data[t1:t2])

        # display computed measures
        for measure in self.heartpy_params['m'].keys():
            print('%s: %f' % (measure, self.heartpy_params['m'][measure]))

    def find_p_peaks(self):
        # Prepare P windows array
        p_windows = np.tile(np.array(np.ones(shape=self.saecg_window_size)), (len(self.beats), 1))

        for i in range(len(self.beats)):
            beat_len = len(self.beats[i])
            self.p_window_center_s[
                i] = beat_len - self.p_window_center * self.fs  # p window center in samples (use as index)
            p_windows[i] = self.beats[i][round(self.p_window_center_s[i] - 0.5 * self.saecg_window_size):round(
                self.p_window_center_s[i] + 0.5 * self.saecg_window_size)]
            #    print('p windows.shape : ', p_windows.shape)
            #    print('p windows : ', p_windows)
            # temp = np.where(p_windows[i] == np.amax(p_windows[i]))[0]
            temp = sp.signal.find_peaks(p_windows[i], distance=len(p_windows[i]))
            if len(temp) > 1:  # If there are equal maxima with same value, take first one
                temp = temp[0]
            #    print('self.p_max_locs : ', self.p_max_locs)
            self.p_max_locs[i] = temp

            #   print('p_max_locs[i] = ', self.p_max_locs[i])
            self.p_max_locs[i] = self.p_max_locs[i] + self.p_window_center_s[i] - round(0.5 * self.saecg_window_size)

        return self
        # for i in range(len(beats)):
        #    print('p peak location = ', self.p_max_locs[i])

    def set_beats(self, beats):
        self.beats = beats

    def get_beats(self):
        return self.beats

    def get_p_max_locs(self):
        return self.p_max_locs

    def set_p_max_locs(self, p_max_locs):
        self.p_max_locs = p_max_locs

    def compute_saecg(self):
        for i in range(len(self.p_max_locs)):
            self.saecg_p_windows[i] = self.get_beats()[i][int(self.p_max_locs[i] - 0.5 * self.saecg_window_size):
                                                          int(self.p_max_locs[i] + 0.5 * self.saecg_window_size)]

        self.saecg_p = np.sum(self.saecg_p_windows, axis=0)
        self.saecg_p = np.true_divide(self.saecg_p, len(self.get_beats()))
        return self

    # Compute the xcorr of SAECG with the array of beats using le and re of P template (specified by user)
    def compute_saecg_xcorr(self, le, re):
        saecg_local = self.saecg_p[int(le * self.fs):int(re * self.fs)]  # saecg cropped with le and re
        for i in range(len(self.get_beats())):
            # percent change data
            pc_change_saecg_local = 100 * np.true_divide(np.diff(saecg_local), saecg_local[1:])
            pc_change_beats_i = 100 * np.true_divide(np.diff(self.beats[i]), self.beats[i][1:])

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
            # # plt.plot(xcorr_p_window)
            # plt.plot(xcorr_p_window[0:self.xcorr_params['max_loc'][i]])
            # plt.plot(self.xcorr_params['le'][i], xcorr_p_window[self.xcorr_params['le'][i]], 'x')
            # plt.plot(self.xcorr_params['max_loc'][i], xcorr_p_window[self.xcorr_params['max_loc'][i]], 'o')
            # plt.show()

            p_window_le = round(0.5 * len(xcorr))  # loc of p window le in samples from start of beat
            print('DM - p_window_le', p_window_le)
            print('DM - 0.5*len(saecg_local)', 0.5 * len(saecg_local))
            print('DM - self.p_window_center_s[i]', self.p_window_center_s[i])
            print('DM - p_window_max_loc=', self.xcorr_params['max_loc'][i])
            print('len(xcorr)', len(xcorr))
            print('len(self.beats[i])', len(self.beats[i]))

            # Set the loc params to be in terms of beat rather than p window
            self.xcorr_params['max_loc'][i] = self.xcorr_params['max_loc'][i] + p_window_le
            self.xcorr_params['le'][i] = self.xcorr_params['le'][i] + p_window_le
            self.xcorr_params['re'][i] = self.xcorr_params['re'][i] + p_window_le

            # print('temp.tolist().index(np.amax(temp)) = ', temp.tolist().index(np.amax(temp)))

        print('DM - self.self.xcorr_params[max_loc][i] : ', self.xcorr_params['max_loc'][i])
        print('DM - self.xcorr_params[max][i] : ', self.xcorr_params['max_loc'][i])
        return self

    def get_saecg_p(self):
        return self.saecg_p

    class IO:
        def __init__(self):
            print('IO - init')

        # Save the passed dict as a csv file
        def dict_as_csv(self, d):
            pd.DataFrame(d).to_csv('test.csv')

        def load_data(self, fs):
            print('IO - load data')
            signal_duration = 240  # signal duration in seconds

            raw_data = np.loadtxt(open('Data/e0103.csv', "rb"), skiprows=1)
            return raw_data[0:signal_duration * fs, 1]


class UIManager:
    def __init__(self, index, dm):
        self.current_index = index
        self.dm = dm
        self.raw_data = dm.raw_data
        self.p_max_locs = dm.p_max_locs
        self.saecg_p = dm.saecg_p
        self.fs = dm.fs
        self.mode = 1  # view mode for beat viewer

        self.fig, self.ax = plt.subplots(nrows=5, ncols=1, gridspec_kw={'height_ratios': [1, 1, 0.1, 0.1, 1]})
        self.fig.tight_layout()  # Configure the layout of UI elements

        self.axprev = plt.axes([0.7, 0.01, 0.08, 0.03])
        self.axnext = plt.axes([0.81, 0.01, 0.08, 0.03])
        self.axanalyze = plt.axes([0.05, 0.01, 0.08, 0.03])
        self.axmode = plt.axes([0.15, 0.01, 0.08, 0.03])

        self.bnext = Button(self.axnext, 'Next')
        self.bnext.on_clicked(lambda x: self.next_button_pushed(x))

        self.bprev = Button(self.axprev, 'Previous')
        self.bprev.on_clicked(lambda x: self.prev_button_pushed(x))

        self.banalyze = Button(self.axanalyze, 'Analyze')
        self.banalyze.on_clicked(lambda x: self.analyze_button_pushed(x))

        self.bmode = Button(self.axmode, 'Mode')
        self.bmode.on_clicked(lambda x: self.mode_button_pushed(x))

        plt.axes(self.ax[0])
        self.ax[0].set_title('Raw Data')
        plt.plot(np.true_divide(range(len(self.raw_data)), self.fs), self.raw_data)
        plt.scatter(np.true_divide(self.dm.heartpy_params['wd']['peaklist'], self.fs),
                    self.raw_data[self.dm.heartpy_params['wd']['peaklist']], c='g')

        plt.axes(self.ax[1])
        title_string = 'Signal Averaged P Wave'
        self.ax[1].set_title(title_string)
        plt.plot(np.true_divide(range(len(self.saecg_p)), self.fs), self.saecg_p)

        # Make the slider for left edge
        plt.axes(self.ax[2])
        self.slider_start = plt.Slider(self.ax[2], 'Start', 0, len(self.saecg_p) / self.fs,
                                       valinit=0.1 * len(self.saecg_p) / self.fs, color='blue')
        self.slider_start.on_changed(lambda x: self.slider_updated(x))
        # Make the slider for right edge
        plt.axes(self.ax[3])
        self.slider_end = plt.Slider(self.ax[3], 'End', 0, len(self.saecg_p) / self.fs,
                                     valinit=0.9 * len(self.saecg_p) / self.fs, color='blue')
        self.slider_end.on_changed(lambda x: self.slider_updated(x))

        plt.axes(self.ax[1])
        self.ax[1].vlines(self.slider_start.val, self.saecg_p[int(self.slider_start.val * self.fs)] - 0.02,
                          self.saecg_p[int(self.slider_start.val * self.fs)] + 0.02, color='red')
        self.ax[1].vlines(self.slider_end.val, self.saecg_p[int(self.slider_end.val * self.fs)] - 0.02,
                          self.saecg_p[int(self.slider_end.val * self.fs)] + 0.02, color='red')

        plt.axes(self.ax[4])
        self.ax[4].clear()
        title_string = 'Beat %d out of %d beats' % (self.get_current_index() + 1, len(dm.get_beats()))
        self.ax[4].set_title(title_string)

        plt.plot(np.true_divide(range(len(self.dm.get_beats()[self.current_index])), self.fs),
                 self.dm.get_beats()[self.current_index])

        plt.plot(self.p_max_locs[self.get_current_index()] / self.fs,
                 self.dm.get_beats()[self.get_current_index()][round(self.p_max_locs[self.get_current_index()])],
                 marker='x')

        plt.plot(self.dm.xcorr_params['max_loc'][self.get_current_index()] / self.fs,
                 self.dm.get_beats()[self.current_index][self.dm.xcorr_params['max_loc'][self.get_current_index()]],
                 'o')

        mng = plt.get_current_fig_manager()
        mng.window.state("zoomed")
        plt.show()

    def refresh_beat_viewer(self):
        if self.get_current_index() >= 0:
            plt.axes(self.ax[4])
            self.ax[4].clear()
            title_string = 'Beat %d out of %d beats' % (self.get_current_index() + 1, len(dm.get_beats()))
            self.ax[4].set_title(title_string)

            if self.mode == 1:
                plt.plot(np.true_divide(range(len(self.dm.get_beats()[self.current_index])), self.fs),
                         self.dm.get_beats()[self.current_index])

                plt.plot(self.p_max_locs[self.get_current_index()] / self.fs,
                         self.dm.get_beats()[self.get_current_index()][
                             round(self.p_max_locs[self.get_current_index()])],
                         marker='x')

                plt.plot(self.dm.xcorr_params['max_loc'][self.get_current_index()] / self.fs,
                         self.dm.get_beats()[self.current_index][
                             self.dm.xcorr_params['max_loc'][self.get_current_index()]],
                         marker='o')

                plt.plot(self.dm.xcorr_params['re'][self.get_current_index()] / self.fs,
                         self.dm.get_beats()[self.current_index][
                             self.dm.xcorr_params['re'][self.get_current_index()]],
                         marker='*')

            elif self.mode == 2:
                # Create the array which is the 3 beats index-1 through index + 1
                display_data = np.append(
                    np.append(self.dm.get_beats()[self.current_index - 1], self.dm.get_beats()[self.current_index]),
                    self.dm.get_beats()[self.current_index + 1])
                plt.plot(np.true_divide(range(len(display_data)), self.fs), display_data)

                plt.plot((len(self.dm.get_beats()[self.current_index - 1]) + self.dm.xcorr_params['max_loc'][
                    self.get_current_index()]) / self.fs,
                         display_data[
                             (len(self.dm.get_beats()[self.current_index - 1]) + self.dm.xcorr_params['max_loc'][
                                 self.get_current_index()])],
                         marker='o')

                plt.plot(
                    (len(self.dm.get_beats()[self.current_index - 1]) + self.p_max_locs[self.current_index]) / self.fs,
                    display_data[
                        len(self.dm.get_beats()[self.current_index - 1]) + self.p_max_locs[self.current_index]],
                    marker='x')

                plt.plot((len(self.dm.get_beats()[self.current_index - 1]) + self.dm.xcorr_params['re'][
                    self.get_current_index()]) / self.fs,
                         self.dm.get_beats()[self.current_index][
                             self.dm.xcorr_params['re'][self.get_current_index()]],
                         marker='*')

            plt.show()

    def next_button_pushed(self, event):
        # print('data type of index : ', type(self.get_current_index()))
        print('self.dm.xcorr_params[max_loc] = ', self.dm.xcorr_params['max_loc'])
        if self.get_current_index() < len(self.dm.get_beats()) - 1:
            plt.axes(self.ax[4])
            self.set_current_index(self.get_current_index() + 1)
            self.ax[4].clear()
            title_string = 'Beat %d out of %d beats' % (self.get_current_index() + 1, len(self.dm.get_beats()))
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
            self.set_current_index(self.get_current_index() - 1)
            self.ax[4].clear()
            title_string = 'Beat %d out of %d beats' % (self.get_current_index() + 1, len(dm.get_beats()))
            self.ax[4].set_title(title_string)
            self.refresh_beat_viewer()

    def analyze_button_pushed(self, event):
        self.dm = self.dm.find_p_peaks()
        self.dm = self.dm.compute_saecg()
        self.dm = self.dm.compute_saecg_xcorr(self.slider_start.val, self.slider_end.val)
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

    def slider_updated(self, event):
        plt.axes(self.ax[1])
        self.ax[1].clear()
        title_string = 'Signal Averaged P Wave'
        self.ax[1].set_title(title_string)
        plt.plot(np.true_divide(range(len(self.saecg_p)), self.fs), self.saecg_p)
        self.ax[1].vlines(self.slider_start.val, self.saecg_p[int(self.slider_start.val * self.fs)] - 0.02,
                          self.saecg_p[int(self.slider_start.val * self.fs)] + 0.02, color='red')
        self.ax[1].vlines(self.slider_end.val, self.saecg_p[int(self.slider_end.val * self.fs)] - 0.02,
                          self.saecg_p[int(self.slider_end.val * self.fs)] + 0.02, color='red')

    def get_current_index(self):
        return self.current_index

    def set_current_index(self, index):
        self.current_index = index


# START HERE
dm = dataManager()
ui = UIManager(1, dm)

plt.show()
