function sum_stat = get_stat(activation, baseline, mask)

    % use ft_freqstatistics to compute t-vals
    cfg = [];
    cfg.channel          = {'all', '-E126', '-E127'};
    cfg.latency          = [-.3 -.05];
    cfg.frequency        = [-5 5]; % relative to terminal f0
    cfg.avgoverfreq      = 'yes';
    cfg.method           = 'analytic';
    cfg.statistic        = 'ft_statfun_actvsblT'; 
    cfg.correctm         = 'no';
    N = size(activation.powspctrm, 1);
    cfg.design  = [ones(1,N) 2*ones(1,N); 1:N, 1:N];
    cfg.ivar = 1; 
    cfg.uvar = 2;
    [stat] = ft_freqstatistics(cfg, activation, baseline);
    
    % mask the t-val array to just specified cluster
    t_vals = stat.stat;
    t_vals(~mask) = 0; % apply mask
    % sum over cluster to get maxsum-like test statistic
    cluster_sum = sum(t_vals(:));
    % and return
    sum_stat = cluster_sum;

end