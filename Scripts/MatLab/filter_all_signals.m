% run signal_filtering or bandpass_filter on all Hb_tot data 
% (see function.m files for input args)
path = "Data/signal_segments/Hb_tot/";
Hbtot_files = dir(fullfile(path, "*.txt"));
Hbtot_filt = {};

for i = 1:numel(Hbtot_files)
    data = fullfile(path, Hbtot_files(i).name);
    Hbtot_filt{i} = bandpass_filter(data, 0.02, 120, 0.001, 1);
end

% run signal_filtering or bandpass_filter on all Hb_oxy data
path = "Data/signal_segments/Hb_oxy/";
Hboxy_files = dir(fullfile(path, "*.txt"));
Hboxy_filt = {};

for i = 1:numel(Hboxy_files)
    data = fullfile(path, Hboxy_files(i).name);
    Hboxy_filt{i} = bandpass_filter(data, 0.02, 120, 0.001, 1);
end

% save cells as files and export to R for analysis



