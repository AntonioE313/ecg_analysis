import heartpy as hp
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
from matplotlib.widgets import Button


class dataManager:
    def __init__(self, beats):
        self.beats = beats
        self.fs = 250
        self.saecg_window_size = int(0.16 * self.fs)  # Watch for rounding error here, need to convert to int
        self.p_max_locs = np.array([[0] * len(beats)])[0]  ##prealocate size beats[i] -> P_max_locs[i]
        self.saecg_xcorr_vals = np.array([[0] * len(beats)])[0]
        self.saecg_p_windows = np.tile(np.array(np.ones(shape=(self.saecg_window_size))), (len(self.beats), 1))
        self.saecg_p = np.array(np.ones(shape=self.saecg_window_size))
        self.find_p_peaks()
        self.compute_saecg()
        self.compute_saecg_xcorr()

    def find_p_peaks(self):
        # Prepare P windows array
        p_window_center = 0.2  # P window centre as in s back from end of beat
        p_windows = np.tile(np.array(np.ones(shape=(self.saecg_window_size))), (len(self.beats), 1))

        for i in range(len(beats)):
            beat_len = len(beats[i])
            p_window_center_s = beat_len - p_window_center * self.fs  # p window center in samples (use as index)
            p_windows[i] = beats[i][round(p_window_center_s - 0.5 * self.saecg_window_size):round(
                p_window_center_s + 0.5 * self.saecg_window_size)]
            #    print('p windows.shape : ', p_windows.shape)
            #     print('p windows : ', p_windows)
            temp = np.where(p_windows[i] == np.amax(p_windows[i]))[0]
            if len(temp) > 1:  # If there are equal maxima with same value, take first one
                temp = temp[0]
            #    print('self.p_max_locs : ', self.p_max_locs)
            self.p_max_locs[i] = temp

            #   print('p_max_locs[i] = ', self.p_max_locs[i])
            self.p_max_locs[i] = self.p_max_locs[i] + p_window_center_s - round(0.5 * self.saecg_window_size)

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
            #        print('i value : ', i)
            #        print(
            # 'self.get_beats()[i][int(self.p_max_locs[i] - 0.5*self.saecg_window_size):
            # int(self.p_max_locs[i] + 0.5*self.saecg_window_size)].shape : ',
            #            self.get_beats()[i][int(self.p_max_locs[i] - 0.5 * self.saecg_window_size):int(
            #                self.p_max_locs[i] + 0.5 * self.saecg_window_size)].shape)
            #        print('self.saecg_window_size : ', self.saecg_window_size)
            #        print('int(self.p_max_locs[i] - 0.5*self.saecg_window_size) : ',
            #              int(self.p_max_locs[i] - 0.5 * self.saecg_window_size))
            #        print('int(self.p_max_locs[i] + 0.5*self.saecg_window_size) : ',
            #              int(self.p_max_locs[i] + 0.5 * self.saecg_window_size))
            #        print('len(self.get_beats()[i]) : ', len(self.get_beats()[i]))

            self.saecg_p_windows[i] = self.get_beats()[i][int(self.p_max_locs[i] - 0.5 * self.saecg_window_size):
                                                          int(self.p_max_locs[i] + 0.5 * self.saecg_window_size)]

        #    print('self.saecg_p_windows[i].dtype : ', self.saecg_p_windows[i].dtype)

        #    print('self.saecg_p_windows : ', self.saecg_p_windows)
        #    print('type(self.saecg_p_windows) : ', type(self.saecg_p_windows))
        #    print('self.saecg_p_windows.shape : ', self.saecg_p_windows.shape)
        self.saecg_p = np.sum(self.saecg_p_windows, axis=0)
        self.saecg_p = np.true_divide(self.saecg_p, len(self.get_beats()))

    def compute_saecg_xcorr(self):
        for i in range(len(self.beats)):
            self.saecg_xcorr_vals[i] = np.amax(sp.signal.correlate(self.saecg_p, self.beats[i]))
        self.saecg_xcorr_vals = np.true_divide(self.saecg_xcorr_vals, len(self.saecg_p))
        #    print('self.saecg_xcorr_vals[i] : ', self.saecg_xcorr_vals[i])
        #   plt.figure()
        #  plt.plot(sp.signal.correlate(self.saecg_p, self.beats[i]))

    def get_saecg_p(self):
        return self.saecg_p


class UIManager:
    def __init__(self, index, beats, p_max_locs, wd, m, raw_data, saecg_p, saecg_xcorr_vals, fs):
        self.current_index = index
        self.beats = beats
        self.wd = wd
        self.m = m
        self.raw_data = raw_data
        self.p_max_locs = p_max_locs
        self.saecg_p = saecg_p
        self.fs = fs
        self.saecg_xcorr_vals = saecg_xcorr_vals
        self.fig, self.ax = plt.subplots(nrows=3, ncols=1)

        self.axprev = plt.axes([0.7, 0.02, 0.08, 0.05])
        self.axnext = plt.axes([0.81, 0.02, 0.08, 0.05])
        self.bnext = Button(self.axnext, 'Next')
        self.bnext.on_clicked(lambda x: self.next_button_pushed(x))
        self.bprev = Button(self.axprev, 'Previous')
        self.bprev.on_clicked(lambda x: self.prev_button_pushed(x))

        plt.axes(self.ax[0])
        self.ax[0].set_title('Raw Data')
        plt.plot(raw_data)
        plt.scatter(wd['peaklist'], raw_data[wd['peaklist']], c='g')
        plt.plot(wd['peaklist'][1:-1], 6 * self.saecg_xcorr_vals, marker='_')

        plt.axes(self.ax[1])
        plt.plot(self.saecg_p)

        plt.axes(self.ax[2])
        self.ax[2].clear()
        title_string = 'Beat %d out of %d beats' % (self.get_current_index() + 1, len(self.get_beats()))
        self.ax[2].set_title(title_string)
        plt.plot(np.true_divide(range(len(self.get_beats()[self.get_current_index()])), self.fs),
                 self.get_beats()[self.get_current_index()])

        mng = plt.get_current_fig_manager()
        mng.window.state("zoomed")
        plt.show()

    def next_button_pushed(self, event):
        print('data type of index : ', type(self.get_current_index()))
        if self.get_current_index() < len(self.get_beats()) - 1:
            self.set_current_index(self.get_current_index() + 1)
            self.ax[2].clear()
            title_string = 'Beat %d out of %d beats' % (self.get_current_index() + 1, len(self.get_beats()))
            self.ax[2].set_title(title_string)
            plt.plot(np.true_divide(range(len(self.get_beats()[self.get_current_index()])), self.fs),
                     self.get_beats()[self.get_current_index()])

            print('Type of self.p_max_locs[self.get_current_index()]: ',
                  type(self.p_max_locs[self.get_current_index()]))
            print('Type of self.get_current_index() : ', type(self.get_current_index()))
            print('data type of index : ', type(self.get_current_index()))
            plt.plot(self.p_max_locs[self.get_current_index()] / self.fs,
                     self.get_beats()[self.get_current_index()][round(self.p_max_locs[self.get_current_index()])],
                     marker='x')
            plt.show()

    def prev_button_pushed(self, event):
        if self.get_current_index() > 0:
            self.set_current_index(self.get_current_index() - 1)
            self.ax[2].clear()
            title_string = 'Beat %d out of %d beats' % (self.get_current_index() + 1, len(self.get_beats()))
            self.ax[2].set_title(title_string)
            plt.plot(np.true_divide(range(len(self.get_beats()[self.get_current_index()])), self.fs),
                     self.get_beats()[self.get_current_index()])
            plt.plot(self.p_max_locs[self.get_current_index()] / self.fs,
                     self.get_beats()[self.get_current_index()][round(self.p_max_locs[self.get_current_index()])],
                     marker='x')
            plt.show()

    def get_current_index(self):
        return self.current_index

    def set_current_index(self, index):
        self.current_index = index

    def set_beats(self, beats):
        self.beats = beats

    def get_beats(self):
        return self.beats


def startup():
    signal_duration = 240  # signal duration in seconds

    sample_rate = 250

    data = np.loadtxt(open('e0103.csv', "rb"), skiprows=1)
    data = data[0:signal_duration * sample_rate, 1]

    plt.figure(figsize=(12, 4))
    plt.plot(data)
    plt.show()

    # run analysis
    wd, m = hp.process(data, sample_rate)
    print(wd['peaklist'])
    print(wd['peaklist'][1])

    print(type(wd['peaklist']))
    peaklist = wd['peaklist']
    beats = []
    for count, value in enumerate(peaklist, start=2):
        print(value, count)
        if count < len(peaklist):
            t1 = peaklist[count - 2]
            t2 = peaklist[count - 1]
            beats.append(data[t1:t2])

    # display computed measures
    for measure in m.keys():
        print('%s: %f' % (measure, m[measure]))

    #  p_max_locs = max(beats)

    return beats, wd, m, data


# START HERE
beats, wd, m, raw_data = startup()
dm = dataManager(beats)
print("The length of beats is: ", len(dm.get_beats()))
ui = UIManager(0, dm.get_beats(), dm.get_p_max_locs(), wd, m, raw_data, dm.get_saecg_p(), dm.saecg_xcorr_vals, dm.fs)

plt.show()