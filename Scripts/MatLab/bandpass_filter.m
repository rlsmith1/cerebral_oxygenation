function bp = bandpass_filter(data, mean_length, cut, low_thresh, high_thresh)

    %% read in data
    data = importdata(data);
    variable = data.data(:,2);
    
    %% average and cut
    Fs = 50; % Hz
    windowSize = mean_length*Fs; % window size needed to get s second rolling average
    a = 1;
    b = (1/windowSize)*ones(1,windowSize);
    variable_movmean = filter(b, a, variable); % mean_length (s) moving average
    variable_movmean_cut = variable_movmean((cut/mean_length):(numel(variable_movmean) - (cut/mean_length))); % take s seconds off each side

    %% filter
    bp = bpfilt(variable_movmean_cut, low_thresh, high_thresh, Fs, false);
    


