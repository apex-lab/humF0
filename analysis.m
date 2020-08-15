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
idx = trialcounts >= 30;
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
cfg.neighbourdist = 3;
neighbors  = ft_prepare_neighbours(cfg);

%% compare baseline period to activation period with permutation test
cfg = [];
cfg.channel          = {'all', '-E126', '-E127'}; % all except eyes
cfg.latency          = [-.3 -.05];
cfg.frequency        = [-5 5]; % relative to terminal f0
cfg.avgoverfreq      = 'yes';
cfg.method           = 'montecarlo';
cfg.statistic        = 'ft_statfun_actvsblT'; % activation vs baseline stat
cfg.correctm         = 'cluster';
cfg.clusteralpha     = 0.025;
cfg.clusterstatistic = 'maxsum';
cfg.alpha            = 0.005;
cfg.numrandomization = 5000;
cfg.neighbours = neighbors;

% specify study design
N = length(allSubj_TFRs);
cfg.design  = [ones(1,N) 2*ones(1,N); 1:N, 1:N];
cfg.ivar = 1; % the 1st row in cfg.design contains IV
cfg.uvar = 2; % and the second contains the subject number (unit variable)

[stat] = ft_freqstatistics(cfg, activationAll, baselineAll);
save('stat.mat', 'stat');

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


%% single plot (without mask)
% baseline the TFR
cfg.baseline = [-.55 -.3];
cfg.baselinetype = 'relative'; % 'relative' computes ERSP
cfg.parameter = 'powspctrm';
[TFR_ersp] = ft_freqbaseline(cfg, TFR);

cfg = [];
cfg.maskstyle = 'saturation';
cfg.ylim = [-40 40];
cfg.xlim = [-.3 -.05];
%cfg.zlim = [.8 ];
cfg.channel      = 'E71';
cfg.elec = ELEC;
figure;
ft_singleplotTFR(cfg, TFR_ersp);  


%% cluster plot (masked)
pos_cluster_pvals = [stat.posclusters(:).prob];
pos_clust = find(pos_cluster_pvals < 0.0025);
pos       = ismember(stat.posclusterslabelmat, pos_clust);

% baseline the TFR
cfg.baseline = [-.55 -.3];
cfg.baselinetype = 'relative'; % 'relative' computes ERSP
cfg.parameter = 'powspctrm';
[TFR_ersp] = ft_freqbaseline(cfg, TFR);
% and then remove everything but activation window
cfg= [];
cfg.channel = {'all', '-E126', '-E127'}; % rm eye electrodes
cfg.latency = [-.3 -.05];
cfg.frequency = [-5 5];
cfg.avgoverfreq = 'yes';
activation_ersp = ft_selectdata(cfg, TFR_ersp);
activation_ersp.powspctrm = activation_ersp.powspctrm .* (stat.posclusterslabelmat == 1);


% now make the plot
timestep      = 0.025; % between time windows for each subplot
sample_count  = length(activation_ersp.time);
j = [-.3:timestep:-.05]; 


[i1,i2] = match_str(activation_ersp.label, stat.label);

for k = 1:(length(j)-1)
   subplot(2,5,k);
   cfg = [];
   cfg.xlim = [j(k) j(k+1)];   % time interval of the subplot
   cfg.zlim = [0 1.2];
   % If a channel is in a to-be-plotted cluster, then
   % the element of pos_int with an index equal to that channel
   % number will be set to 1 (otherwise 0).

   % Next, check which channels are in the clusters over the
   % entire time interval of interest.
   pos_int = zeros(numel(activation_ersp.label),1);

   cfg.highlight   = 'on';
   % Get the index of the to-be-highlighted channel
   cfg.highlightchannel = 'E71';
   cfg.highlightcolor = [1 0 0];
   cfg.highlightsymbol = 'o';
   cfg.highlightsize = 1;
   cfg.markersize = 1;
   cfg.comment     = 'xlim';
   cfg.commentpos  = 'title';
   cfg.elec      = ELEC;
   cfg.interactive = 'no';
   cfg.interplimits = 'head';
   ft_topoplotER(cfg, activation_ersp);
end






















