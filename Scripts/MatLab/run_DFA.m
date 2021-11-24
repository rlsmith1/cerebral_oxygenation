
%% run DFA on all filtered signals

% Hb_tot
N = numel(Hbtot_filt);
Hbtot_alphas = zeros(N, 1);
for i = 1:N
    signal = Hbtot_filt{i}(180:(numel(Hbtot_filt{i})-180),:); % chop three minutes off beg & end
    [D, alpha] = DFA_main(signal);
    Hbtot_alphas(i,1) = alpha;
end

% Hb_oxy
N = numel(Hboxy_filt);
Hboxy_alphas = zeros(N, 1);
for i = 1:N
    signal = Hboxy_filt{i}(180:(numel(Hboxy_filt{i})-180),:); % chop three minutes off beg & end
    [D, alpha] = DFA_main(signal);
    Hboxy_alphas(i,1) = alpha;
end


