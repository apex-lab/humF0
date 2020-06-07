FT_DIRECTORY = '~/repos/fieldtrip';
addpath(FT_DIRECTORY)

F0 = [166.1224 % the average terminal f0 of each subject, in order
  105.6851
  103.4111
  108.3090
  270.9184
   75.0000
  181.8367
  113.5569
  124.0525
   98.4694
  124.3163
   97.2303
  200.0902
  199.3586
  194.8688
];

SUBJECTS = 7:21; % subject numbers
ELEC = strcat(FT_DIRECTORY, '/template/electrode/GSN-HydroCel-129.sfp');

%% compute time-frequency representations (TFRs) for each subject
allSubj_TFRs = cell(1, length(SUBJECTS));
allSubj_baseline = cell(1, length(SUBJECTS));
allSubj_activation = cell(1, length(SUBJECTS));
trialcounts = 1:length(SUBJECTS); % placeholder

for i = 1:length(SUBJECTS)
    try 
        data = preprocess(SUBJECTS(i));
        trialcounts(i) = length(data.trial);
        [TFR, TFR_baseline, TFR_activation] = compute_TFRs(data, F0(i));
        allSubj_TFRs{i} = TFR;
        allSubj_baseline{i} = TFR_baseline;
        allSubj_activation{i} = TFR_activation;
    catch % error in which there are no trials after artifact rejection
        trialcounts(i) = 0;
    end
    clear data TFR TFR_baseline TFR_activation
end

save('tfr.mat', 'allSubj_TFRs', 'allSubj_baseline', 'allSubj_activation');
save('trial_counts.mat', 'trialcounts');
clear allSubj_TFRs allSubj_baseline allSubj_activation

%% aggregate across subjects
load('tfr.mat', 'allSubj_TFRs', 'allSubj_baseline', 'allSubj_activation');
load('trial_counts.mat', 'trialcounts');
% remove subjects with no good trials
idx = trialcounts >= 15;
allSubj_TFRs = allSubj_TFRs(idx);
allSubj_baseline = allSubj_baseline(idx);
allSubj_activation = allSubj_activation(idx);
% and aggregate
cfg = [];
cfg.keepindividual = 'no';
TFR = ft_freqgrandaverage(cfg, allSubj_TFRs{:});
cfg.keepindividual = 'yes';
baselineAll = ft_freqgrandaverage(cfg, allSubj_baseline{:});
activationAll = ft_freqgrandaverage(cfg, allSubj_activation{:});

%% compute neighbors from average electrode placement
cfg         = [];
cfg.elec    = ELEC;
cfg.method  = 'distance';
neighbors  = ft_prepare_neighbours(cfg);

%% compare baseline period to activation period with permutation test
cfg = [];
cfg.channel          = {'all', '-E126', '-E127'}; % all except eyes
cfg.latency          = 'all';
cfg.frequency        = [-40 40]; % relative to terminal f0
cfg.method           = 'montecarlo';
cfg.statistic        = 'ft_statfun_actvsblT'; % activation vs baseline stat
cfg.correctm         = 'cluster';
cfg.clusteralpha     = 0.05;
cfg.clusterstatistic = 'maxsum';
cfg.alpha            = 0.005;
cfg.numrandomization = 'all';
cfg.neighbours = neighbors;

% specify study design
N = length(allSubj_TFRs);
cfg.design  = [ones(1,N) 2*ones(1,N); 1:N, 1:N];
cfg.ivar = 1; % the 1st row in cfg.design contains IV
cfg.uvar = 2; % and the second contains the subject number (unit variable)

[stat] = ft_freqstatistics(cfg, activationAll, baselineAll);


%% plot group average TFR without mask
cfg = [];
cfg.maskstyle    = 'saturation';
%cfg.zlim           = [1 5];
cfg.baseline       = [-.55 -.3];
cfg.ylim           = [-40 40];
cfg.xlim           = [-.3 -0.05];
cfg.baselinetype   = 'relchange';
cfg.colormap = jet(30);
cfg.channel      = {'all', '-E126', '-E127'};
cfg.elec = ELEC;
figure;
ft_multiplotTFR(cfg,TFR);

%% 
cfg = [];
cfg.maskstyle    = 'saturation';
%cfg.zlim           = [1 5];
cfg.baseline       = [-.55 -.3];
cfg.ylim           = [-40 40];
cfg.xlim           = [-.3 -0.05];
cfg.baselinetype   = 'relative';
cfg.colormap = jet(30);
cfg.channel      = 'Cz';
cfg.elec = ELEC;
figure;
ft_singleplotTFR(cfg,TFR);  

%% plot clusters
cfg = [];
cfg.alpha     = 0.05;
cfg.parameter = 'stat';
%cfg.zlim      = [-4 4];
cfg.elec    = ELEC;
ft_clusterplot(cfg, stat);

%% plot with mask
% baseline the TFR
cfg.baseline = [-.55 -.3];
cfg.baselinetype = 'relative'; % 'relative' computes ERSP
[TFR_ersp] = ft_freqbaseline(cfg, TFR);
% and then remove everything but activation window
cfg= [];
cfg.channel = {'all', '-E126', '-E127'}; % rm eye electrodes
cfg.latency = [-.3 -.05];
cfg.frequency = [-40 40];
activation_ersp = ft_selectdata(cfg, TFR_ersp);

% apply mask, only plotting first significant cluster
% (since we only found one)
activation_ersp.powspctrm = activation_ersp.powspctrm .* (stat.posclusterslabelmat == 1);
% and plot
cfg = [];
cfg.zlim           = [0 3];
cfg.colormap = jet(30);
cfg.channel = 'all';
cfg.elec = ELEC;
figure;
ft_multiplotTFR(cfg,activation_ersp);

%% plot clusters
cfg = [];
cfg.alpha     = 0.01;
cfg.parameter = 'prob';
cfg.colormap = jet(30);
%cfg.zlim      = [-4 4];
cfg.elec    = ELEC;
ft_clusterplot(cfg, stat);









