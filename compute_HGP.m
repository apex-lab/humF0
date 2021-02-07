function activation_ersp = compute_HGP(data)
% Computes log of high gamma power time series

    % downsample data (to save memory)
    cfg = [];
    cfg.resamplefs = 4000; % Hz, new sampling rate
    cfg.feedback = 'no';
    data = ft_resampledata(cfg, data);
    
    % compute high gamma power
    cfg            = [];
    cfg.method     = 'tfr';
    cfg.keeptrials = 'yes';
    cfg.toi        = -0.6:0.002:0; 
    cfg.foi        = [80 85 90 95 100 105 110 115 125 130 135 140 145 ...
        150 155 160 165  ... 
        170 175 185 190]; % 80 - 190 Hz, leave out harmonics of 60 Hz
    cfg.width      = 0.05; % 50 ms width
    cfg.output     = 'pow';
    cfg.keeptrials = 'no';
    TFR = ft_freqanalysis(cfg, data);
    
    % apply baseline
    cfg.baseline = [-.55 -.3];
    cfg.baselinetype = 'relative'; % 'relative' computes ERSP
    cfg.parameter = 'powspctrm';
    TFR_ersp = ft_freqbaseline(cfg, TFR);
    
    % and then remove everything but activation window
    cfg= [];
    cfg.channel = {'all', '-E126', '-E127'}; % rm eye electrodes
    cfg.latency = [-.3 -.05];
    cfg.avgoverfreq = 'yes'; % and average over frequencies
    activation_ersp = ft_selectdata(cfg, TFR_ersp);
    
end