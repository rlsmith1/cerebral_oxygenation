% export Hb_tot signals
for i = 1:numel(Hbtot_filt)
    path = "Data/filtered_signals/0.2s_mean_0.001_1_filt/cerebral_missing_UM/Hb_tot";
    writematrix(Hbtot_filt{i}, fullfile(path, Hbtot_files(i).name));
end

% export Hb_oxy signals
for i = 1:numel(Hboxy_filt)
    path = "Data/filtered_signals/0.2s_mean_0.001_1_filt/cerebral_missing_UM/Hb_oxy";
    writematrix(Hboxy_filt{i}, fullfile(path, Hboxy_files(i).name));
end