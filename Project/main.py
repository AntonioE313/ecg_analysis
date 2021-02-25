import heartpy as hp
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
from matplotlib.widgets import Button
from matplotlib.widgets import Slider


class dataManager:
    def __init__(self, beats, raw_data):
        self.beats = beats
        self.raw_data = raw_data
        self.fs = 250
        self.saecg_window_size = int(0.2 * self.fs)  # Watch for rounding error here, need to convert to int
        self.p_max_locs = np.array([[0] * len(beats)])[0]  # Pre-allocate size beats[i] -> P_max_locs[i]
        self.saecg_xcorr_max_vals = np.array([[0] * len(beats)])[0]
        self.saecg_p_windows = np.tile(np.array(np.ones(shape=(self.saecg_window_size))), (len(self.beats), 1))
        self.saecg_p = np.array(np.ones(shape=self.saecg_window_size))
        self.find_p_peaks()
        self.compute_saecg()

    def find_p_peaks(self):
        # Prepare P windows array
        p_window_center = 0.15  # P window centre as in seconds back from end of beat
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
            self.saecg_p_windows[i] = self.get_beats()[i][int(self.p_max_locs[i] - 0.5 * self.saecg_window_size):
                                                          int(self.p_max_locs[i] + 0.5 * self.saecg_window_size)]

        self.saecg_p = np.sum(self.saecg_p_windows, axis=0)
        self.saecg_p = np.true_divide(self.saecg_p, len(self.get_beats()))
        return self

    def compute_saecg_xcorr(self):
        for i in range(len(self.beats)):
            # percent change data
            temp = np.correlate(np.true_divide(np.diff(self.saecg_p),
                                               self.saecg_p[1:]),
                                np.true_divide(np.diff(self.beats[i]), self.beats[i][1:]), mode='same')
            self.saecg_xcorr_max_vals[i] = np.amax(temp)

            if i < 4:
                print('temp = ', temp)
                plt.figure(i*100 + i*10 + i)
                print(i)
                plt.plot(np.true_divide(range(len(temp)), self.fs), temp)
                plt.show()
    #For some reason saecg_xcorr_max_vals is not keeping value when leaving for loop. Scope issue with multiple copies?
    #I didn't change anything and now it works?
        self.saecg_xcorr_max_vals = np.true_divide(self.saecg_xcorr_max_vals, len(self.saecg_p))
        print('self.saecg_xcorr_max_vals : ', self.saecg_xcorr_max_vals)
        return self

    def get_saecg_p(self):
        return self.saecg_p


class UIManager:
    # ui = UIManager(0, dm.get_beats(), dm.get_p_max_locs(), wd, m, raw_data, dm.get_saecg_p(), dm.saecg_xcorr_vals, dm.fs)
    def __init__(self, index, dm, wd, m):
        self.current_index = index
        self.dm = dm
        self.wd = wd
        self.m = m
        self.raw_data = dm.raw_data
        self.p_max_locs = dm.p_max_locs
        self.saecg_p = dm.saecg_p
        self.fs = dm.fs

        self.fig, self.ax = plt.subplots(nrows=5, ncols=1, gridspec_kw={'height_ratios':[1, 1, 0.1, 0.1, 1]})
        self.fig.tight_layout()  # Configure the layout of UI elements

        self.axprev = plt.axes([0.7, 0.01, 0.08, 0.03])
        self.axnext = plt.axes([0.81, 0.01, 0.08, 0.03])
        self.axanalyze = plt.axes([0.05, 0.01, 0.08, 0.03])

        self.bnext = Button(self.axnext, 'Next')
        self.bnext.on_clicked(lambda x: self.next_button_pushed(x))

        self.bprev = Button(self.axprev, 'Previous')
        self.bprev.on_clicked(lambda x: self.prev_button_pushed(x))

        self.banalyze = Button(self.axanalyze, 'Analyze')
        self.banalyze.on_clicked(lambda x: self.analyze_button_pushed(x))

        plt.axes(self.ax[0])
        self.ax[0].set_title('Raw Data')
        plt.plot(np.true_divide(range(len(raw_data)), self.fs), raw_data)
        plt.scatter(np.true_divide(wd['peaklist'], self.fs), raw_data[wd['peaklist']], c='g')

        plt.axes(self.ax[1])
        title_string = 'Signal Averaged P Wave'
        self.ax[1].set_title(title_string)
        plt.plot(np.true_divide(range(len(self.saecg_p)), self.fs), self.saecg_p)


        # Make the slider for left edge
        plt.axes(self.ax[2])
        self.slider_start = plt.Slider(self.ax[2], 'Start', 0, len(self.saecg_p)/self.fs, valinit=0.1*len(self.saecg_p)/self.fs, color='blue')
        self.slider_start.on_changed(lambda x: self.slider_updated(x))
        # Make the slider for right edge
        plt.axes(self.ax[3])
        self.slider_end = plt.Slider(self.ax[3], 'End', 0, len(self.saecg_p) / self.fs, valinit=0.9*len(self.saecg_p)/self.fs, color='blue')
        self.slider_end.on_changed(lambda x: self.slider_updated(x))

        plt.axes(self.ax[1])
        self.ax[1].vlines(self.slider_start.val, self.saecg_p[int(self.slider_start.val*self.fs)]-0.02,
                        self.saecg_p[int(self.slider_start.val*self.fs)]+0.02, color='red')
        self.ax[1].vlines(self.slider_end.val, self.saecg_p[int(self.slider_end.val*self.fs)] - 0.02,
                          self.saecg_p[int(self.slider_end.val*self.fs)] + 0.02, color='red')

        plt.axes(self.ax[4])
        self.ax[4].clear()
        title_string = 'Beat %d out of %d beats' % (self.get_current_index() + 1, len(dm.get_beats()))
        self.ax[4].set_title(title_string)

        # Create the array which is the 3 beats index-1 through index + 1
        display_data = np.append(np.append(self.dm.get_beats()[self.current_index-1], self.dm.get_beats()[self.current_index]),
                                 self.dm.get_beats()[self.current_index+1])
        plt.plot(np.true_divide(range(len(display_data)), self.fs), display_data)

        plt.plot(self.p_max_locs[self.get_current_index()] / self.fs,
                 self.dm.get_beats()[self.get_current_index()][round(self.p_max_locs[self.get_current_index()])],
                 marker='x')

        mng = plt.get_current_fig_manager()
        mng.window.state("zoomed")
        plt.show()

    def next_button_pushed(self, event):
        print('data type of index : ', type(self.get_current_index()))
        if self.get_current_index() < len(self.dm.get_beats()) - 1:
            plt.axes(self.ax[4])
            self.set_current_index(self.get_current_index() + 1)
            self.ax[4].clear()
            title_string = 'Beat %d out of %d beats' % (self.get_current_index() + 1, len(self.dm.get_beats()))
            self.ax[4].set_title(title_string)

            # Create the array which is the 3 beats index-1 through index + 1
            display_data = np.append(
                np.append(self.dm.get_beats()[self.current_index - 1], self.dm.get_beats()[self.current_index]),
                self.dm.get_beats()[self.current_index + 1])
            plt.plot(np.true_divide(range(len(display_data)), self.fs), display_data)

            print('Type of self.p_max_locs[self.get_current_index()]: ',
                  type(self.p_max_locs[self.get_current_index()]))
            print('Type of self.get_current_index() : ', type(self.get_current_index()))
            print('data type of index : ', type(self.get_current_index()))
            plt.plot(self.p_max_locs[self.get_current_index()] / self.fs,
                     self.dm.get_beats()[self.get_current_index()][round(self.p_max_locs[self.get_current_index()])],
                     marker='x')
            plt.show()

    def prev_button_pushed(self, event):
        if self.get_current_index() > 0:
            plt.axes(self.ax[4])
            self.set_current_index(self.get_current_index() - 1)
            self.ax[4].clear()
            title_string = 'Beat %d out of %d beats' % (self.get_current_index() + 1, len(dm.get_beats()))
            self.ax[4].set_title(title_string)

            # Create the array which is the 3 beats index-1 through index + 1
            display_data = np.append(
                np.append(self.dm.get_beats()[self.current_index - 1], self.dm.get_beats()[self.current_index]),
                self.dm.get_beats()[self.current_index + 1])
            plt.plot(np.true_divide(range(len(display_data)), self.fs), display_data)

            plt.plot(self.p_max_locs[self.get_current_index()] / self.fs,
                     self.dm.get_beats()[self.get_current_index()][round(self.p_max_locs[self.get_current_index()])],
                     marker='x')
            plt.show()

    def analyze_button_pushed(self, event):
        temp = self.saecg_p[int(self.slider_start.val * self.fs):int(self.slider_end.val * self.fs)]
        dm.saecg_p = temp
        self.dm = self.dm.compute_saecg()
        self.dm = self.dm.compute_saecg_xcorr()

    def slider_updated(self, event):
        plt.axes(self.ax[1])
        self.ax[1].clear()
        title_string = 'Signal Averaged P Wave'
        self.ax[1].set_title(title_string)
        plt.plot(np.true_divide(range(len(self.saecg_p)), self.fs), self.saecg_p)
        self.ax[1].vlines(self.slider_start.val, self.saecg_p[int(self.slider_start.val*self.fs)]-0.02,
                        self.saecg_p[int(self.slider_start.val*self.fs)]+0.02, color='red')
        self.ax[1].vlines(self.slider_end.val, self.saecg_p[int(self.slider_end.val*self.fs)] - 0.02,
                          self.saecg_p[int(self.slider_end.val*self.fs)] + 0.02, color='red')

    def get_current_index(self):
        return self.current_index

    def set_current_index(self, index):
        self.current_index = index


def startup():
    signal_duration = 240  # signal duration in seconds

    sample_rate = 250

    data = np.loadtxt(open('Data/e0103.csv', "rb"), skiprows=1)
    data = data[0:signal_duration * sample_rate, 1]

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
dm = dataManager(beats, raw_data)
print("The length of beats is: ", len(dm.get_beats()))
# ui = UIManager(0, dm.get_beats(), dm.get_p_max_locs(), wd, m, raw_data, dm.get_saecg_p(), dm.saecg_xcorr_vals, dm.fs)
ui = UIManager(1, dm, wd, m)

plt.show()
