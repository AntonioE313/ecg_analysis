# Config file for easily changing app behaviour
# data_filepath = 'Data/e0103.csv'
data_filepath = 'Data/paf-prediction-challenge-database-1.0.0/n08'
signal_duration = 400  # Signal duration in seconds

p_window_center = 0.15  # Initial P window centre as in seconds back from peak of QRS

premature_signal_coefficient = 0.99  # Threshold for RR interval i in terms of RR interval (i-1) for premature state=1

upsample_coefficient = 4 # Upsample up this factor for better use interfacing with heartpy