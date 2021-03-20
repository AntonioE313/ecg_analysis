# Config file for easily changing app behaviour
data_filepath = 'Data/e0103.csv'
fs = 250
signal_duration = 60  # Signal duration in seconds

p_window_center = 0.15  # Initial P window centre as in seconds back from peak of QRS
saecg_window_size = int(0.2 * fs)  # Watch for rounding error here converting to int

premature_signal_coefficient = 0.99  # Threshold for RR interval i in terms of RR interval (i-1) for premature state=1
