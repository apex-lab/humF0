function stat = one_permutation(f0s, trialcounts, chans, mask)

    % loops over subjects w/ high enough trial counts
    subs = 1:length(trialcounts);
    idx = trialcounts >= 30;
    subs = subs(idx); 
    perm_activations = cell(1, length(subs));
    perm_baselines = cell(1, length(subs));
    for i = subs
        % load and select chans from preprocessed data
        load(['deriv/sub' num2str(i) 'cleaned.mat'], 'data');
        cfg = [];
        cfg.channel = chans; 
        data = ft_selectdata(cfg, data);
        % compute TFRs centered around random F0
        [~, TFR_baseline, TFR_activation] = compute_TFRs(data, f0s(i));
        perm_baselines{i} = TFR_baseline;
        perm_activations{i} = TFR_activation;
    end
    
    % compute test statistic for specified cluster 
    stat = get_stat(perm_activations, perm_baselines, mask);


end