signal_path = "Data/signal_segments/Hb_tot/TM0001CM01.txt";

mean_length = 0.2;
cut = 120;
low_thresh = 0.001;
high_thresh = 1;

 %% read in data
    data = importdata(signal_path);
    variable = data.data(:,2);
    
    %% average and cut
    Fs = 50; % Hz
    windowSize = mean_length*Fs; % window size needed to get s second rolling average
    a = 1;
    b = (1/windowSize)*ones(1,windowSize);
    variable_movmean = filter(b, a, variable); % mean_length (s) moving average
    variable_movmean_cut = variable_movmean((cut/mean_length):(numel(variable_movmean) - (cut/mean_length))); % take s seconds off each side

    %% filter
    bp = bpfilt(variable_movmean_cut, low_thresh, high_thresh, Fs);