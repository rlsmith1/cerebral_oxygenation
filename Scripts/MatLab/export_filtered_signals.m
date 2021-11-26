% export Hb_tot signals
for i = 1:numel(Hbtot_filt)
    path = "Data/filtered_signals/muscle_0.2s_mean_0.001_1_filt/Hb_tot";
    writematrix(Hbtot_filt{i}, fullfile(path, Hbtot_files(i).name));
end

% export Hb_oxy signals
for i = 1:numel(Hboxy_filt)
    path = "Data/filtered_signals/muscle_0.2s_mean_0.001_1_filt/Hb_oxy";
    writematrix(Hboxy_filt{i}, fullfile(path, Hboxy_files(i).name));
end