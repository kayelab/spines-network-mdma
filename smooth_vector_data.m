function [smoothed_data_vector]=smooth_vector_data(time_series,windowSize)

% windowSize=4 in initial run of pupil smoothing for finding maxima and
% minima
% 
b = (1/windowSize)*ones(1,windowSize);
a = 1;
smoothed_data_vector =filtfilt(b,a,time_series);


