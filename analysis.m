FT_DIRECTORY = '~/repos/fieldtrip';
addpath(FT_DIRECTORY)
ft_defaults;


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
        save(['deriv/sub' num2str(i) 'cleaned.mat'], 'data');
        trialcounts(i) = length(data.trial);
        [TFR, TFR_baseline, TFR_activation] = compute_TFRs(data, F0(i), 60, 1);
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
cfg.clusteralpha     = 0.05/2; % .05 but two tailed
cfg.clusterstatistic = 'maxsum';
cfg.alpha            = 0.0025; % .0025 is just .005 but two-tailed
cfg.numrandomization = 5000;
cfg.neighbours = neighbors;

% specify study design
N = length(allSubj_TFRs);
cfg.design  = [ones(1,N) 2*ones(1,N); 1:N, 1:N];
cfg.ivar = 1; % the 1st row in cfg.design contains IV
cfg.uvar = 2; % and the second contains the subject number (unit variable)

[stat] = ft_freqstatistics(cfg, activationAll, baselineAll);
%save('stat.mat', 'stat');
  
pos_cluster_pvals = [stat.posclusters(:).prob]*2

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
activation_ersp.powspctrm(~pos) = 0;

% now make the plot
timestep      = 0.025; % between time windows for each subplot
sample_count  = length(activation_ersp.time);
j = [-.3:timestep:-.05]; 


[i1,i2] = match_str(activation_ersp.label, stat.label);

f = figure;
[a, ~] = tight_subplot(4, 5, 0, 0, 0);
for k = 1:(length(j)-1)
   set(f, 'currentaxes', a(k));
   cfg = [];
   cfg.xlim = [j(k) j(k+1)];   % time interval of the subplot
   cfg.zlim = [.8 1.5];

   cfg.highlight   = 'on';
   cfg.highlightchannel = 'E71';
   cfg.highlightcolor = [1 0 0];
   cfg.highlightsymbol = 'x';
   cfg.highlightsize = 7;
   cfg.markersize = 2;
   cfg.elec      = ELEC;
   cfg.interactive = 'no';
   cfg.interplimits = 'head';
   cfg.comment = 'no';
   ft_topoplotER(cfg, activation_ersp);
   beg = num2str(j(k)*1000);
   en = num2str(j(k+1)*1000);
   txt = [beg ' to ' en ' ms'];
   title(txt, 'FontSize', 12);
end

% plot (unmasked) TFR of representative electrode below clusterplot 
subplot(4, 5, 11:20);
cfg = [];
cfg.maskstyle = 'saturation';
cfg.ylim = [-40 40];
cfg.xlim = [-.3 -.05];
%cfg.zlim = [.8 1.5];
cfg.channel = 'E71';
cfg.elec = ELEC;
cfg.colorbar = 'yes';
cfg.title = 'Occipitoparietal Electrode E71';
ft_singleplotTFR(cfg, TFR_ersp); 
xlabel('Time (sec) Relative to Voicing Onset', 'FontSize', 15);
ylabel('Frequency (Hz) Relative to F0', 'FontSize', 15);
title('Occipitoparietal Electrode E71', 'FontSize', 15);

%% save figure (adjust dimensions in gui window manually first)
print('f0_cluster_laplacian.png','-dpng','-r500');

%% now run permutation test over F0s
% we'll only use channels in original clusters for tractibility
chan_idx = find(sum(squeeze(pos), 2) > 0);
chans = stat.label(chan_idx);
mask = pos(chan_idx, :, :);

n_perms = 1000;
stats = 1:n_perms;

subs = 1:length(SUBJECTS);
idx = trialcounts >= 30;
subs = subs(idx); 

parpool(5)
parfor perm = 1:n_perms
    
    % randomly shuffle F0s
    f0s = unifrnd(70, 300, [1, length(subs)])
        
    % loops over subjects w/ high enough trial counts
    perm_activations = cell(1, length(subs));
    perm_baselines = cell(1, length(subs));
    for i = 1:length(subs)
        j = subs(i);
        % load and select chans from preprocessed data
        a = load(['deriv/sub' num2str(j) 'cleaned.mat'], 'data');
        data = a.data;
        cfg = [];
        cfg.channel = chans; 
        data = ft_selectdata(cfg, data);
        % compute TFRs centered around random F0
        [~, TFR_baseline, TFR_activation] = compute_TFRs(data, f0s(i), 6, 0);
        perm_baselines{i} = TFR_baseline;
        perm_activations{i} = TFR_activation;
    end
    
    % compute test statistic for specified cluster 
    cfg = [];
    cfg.keepindividual = 'yes';
    perm_baselines = ft_freqgrandaverage(cfg, perm_baselines{:});
    perm_activations = ft_freqgrandaverage(cfg, perm_activations{:});
    stats(perm) = get_stat(perm_activations, perm_baselines, mask);
    
end

% terminate parallel workers
poolobj = gcp('nocreate');
delete(poolobj);

% compare observed to permutation null
obs = stat.posclusters(1).clusterstat;
p_f0 = mean(obs < stats)

%% compute high gamma power
hgp = cell(1, length(subs));
for i = 1:length(subs)
    j = subs(i);
    a = load(['deriv/sub' num2str(j) 'cleaned.mat'], 'data');
    data = a.data;
    hgp{i} = compute_HGP(data);
end

cfg = [];
[grandavg] = ft_freqgrandaverage(cfg, hgp{:});


%% plot high gamma power
timestep      = 0.025; % between time windows for each subplot
sample_count  = length(grandavg.time);
j = [-.3:timestep:-.05]; 
[i1,i2] = match_str(grandavg.label, stat.label);
f = figure;
[a, ~] = tight_subplot(3, 5, 0, 0, 0);
for k = 1:(length(j)-1)
   set(f, 'currentaxes', a(k));
   cfg = [];
   cfg.xlim = [j(k) j(k+1)];   % time interval of the subplot
   cfg.zlim = [.8 1.65];
   cfg.elec      = ELEC;
   cfg.interactive = 'no';
   cfg.interplimits = 'head';
   cfg.comment = 'no';
   ft_topoplotER(cfg, grandavg);
   beg = num2str(j(k)*1000);
   en = num2str(j(k+1)*1000);
   txt = [beg ' to ' en ' ms'];
   title(txt, 'FontSize', 12);
end

% add shared colorbar below topoplots
ax = subplot(3, 5, 11:15);
c = colorbar('north');
caxis([.8 1.65]);
ax.Visible = 'off';

%% estimate stimulus onset asyncrony (SOA)
soas = 1:length(subs);
for i = 1:length(subs)
    j = subs(i);
    a = load(['deriv/sub' num2str(j) 'cleaned.mat'], 'data');
    data = a.data;
    % get onsets from fieldtrip struct
    onsets = data.sampleinfo(:,1) / data.fsample;
    soa = diff(onsets);
    % remove SOAs that are way too big (cause of dropped trials)
    soa = soa(soa < 4);
    % and compute mean SOA for this subject
    soas(i) = mean(soa);
end
mean_soa = mean(soas)
std_soa = std(soas)

















