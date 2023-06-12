function filtered_signal = signal_filtering(data, mean_length, cut, threshold)

    %% read in data
    data = importdata(data);
    variable = data.data(:,2);
    
    %% average and cut
    sampling_freq = 50; % Hz
    % mean_length = 1; % rolling mean length in seconds
    windowSize = mean_length*sampling_freq; % window size needed to get s second rolling average
    a = 1;
    b = (1/windowSize)*ones(1,windowSize);
    variable_movmean = filter(b, a, variable); % mean_length (s) moving average
    
    % cut = 120; % amount of time to take off signal, in seconds
    variable_movmean_cut = variable_movmean((cut/mean_length):(numel(variable_movmean) - (cut/mean_length))); % take s seconds off each side
    % if mod(numel(variable_movmean_cut), 2) == 1
    %     variable_movmean_cut = variable_movmean_cut(2:numel(variable_movmean_cut));
    % end % make sure length is even
    
    %% filter
    % threshold = 0.01; % Hz
    % sampling_rate = mean_length; % after moving average
    % filter_order = 10;
    % cutoff = threshold*sampling_rate*2*pi; % analog filter; must be expressed in rad/s
    % Nyquist = threshold/(sampling_rate/2);
    % [b, a] = butter(filter_order, cutoff, 'high');
    % variable_movmean_cut_filt = filter(b, a, variable_movmean_cut);
    
    %% plot signals
    % figure()
    % hold on
    % plot(variable_movmean_cut)
    % plot(variable_movmean_cut_filt)
    % hold off
    % 
    % figure()
    filtered_signal = highpass(variable_movmean_cut, threshold, sampling_freq);
    
    %% FFT
    % X = variable_movmean_cut_filt; % Hbtot signal
    % Fs = sampling_rate;
    % T = 1/Fs; % sampling period
    % L = numel(X); % length of signal
    % t = (0:L-1)*T; % time vector
    % Y = fft(variable_movmean_cut); % FFT
    % P2 = abs(Y/L);
    % P1 = P2(1:L/2+1);
    % P1(2:end-1) = 2*P1(2:end-1);
    % f = Fs*(0:(L/2))/L;
    % figure()
    % plot(f,P1)


