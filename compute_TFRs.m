function [TFR, TFR_baseline, TFR_activation] = compute_TFRs(data, f0)

%% time-frequency analysis (wavelet)
cfg                  = [];
cfg.channel          = 'all';
cfg.method           = 'wavelet';
wavelet_width        = round(0.05*f0*pi); % wavelet length of ~50ms for F0
cfg.width            = wavelet_width;
cfg.output           = 'pow';
f0                   = round(f0); % so we use an integer frequency
cfg.foi              = (f0 - 60):1:(f0 + 60);
cfg.toi              = -0.6:0.002:0;   
cfg.keeptrials       = 'no';
TFR = ft_freqanalysis(cfg, data); 
% normalize Hz by F0 for later averaging across S's
TFR.freq = -60:1:60; % normalize Hz by F0 for later averaging across S's

%% now compute seperate TFRs on baseline and activation periods
cfg.toi = -0.55:0.002:-0.3;
TFR_baseline = ft_freqanalysis(cfg, data); 
TFR_baseline.freq = -60:1:60;

cfg.toi = -0.3:0.002:-0.05;
TFR_activation = ft_freqanalysis(cfg, data); 
TFR_activation.freq = -60:1:60;

% and set time axis to be equal on both for when we compute stats later
TFR_baseline.time = TFR_activation.time;

end