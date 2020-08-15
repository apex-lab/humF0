function data = preprocess(i) 

%% load data
FT_DIRECTORY = '~/repos/fieldtrip';

s = num2str(i);
while length(s) < 3
    s = strcat('0', s);
end
EEGFILE = strcat('~/Documents/EGG/slowHum/P', s, '.set');
ARTFILE = strcat('~/Documents/EGG/slowHum/P', s, '_art.mat');

ELEC = strcat(FT_DIRECTORY, '/template/electrode/GSN-HydroCel-129.sfp');

addpath(FT_DIRECTORY)
addpath(strcat(FT_DIRECTORY, '/fileio'))

%% define trials for artifact rejection
% sometimes the EGG creates electrical noise as the laryngeal muscles move 
% vertically prior to the first glottis closure, so we exclude the first
% 100 ms before glottis closure from artifact detection. Its important
% then that, when we analyze our data later, that the evidence for our
% effect (e.g. time frequency cluster) does not extend into that window,
% since there will be some leakage into our window of interest depending on
% wavelet length. 
cfg = [];
cfg.datafile                = EEGFILE;
cfg.headerfile              = EEGFILE;
cfg.trialfun                = 'ft_trialfun_general'; % default
cfg.trialdef.eventtype      = 'trigger';
cfg.trialdef.eventvalue     = 'TGCI'; 
cfg.trialdef.prestim        = 0.6; % in seconds
cfg.trialdef.poststim       = -0.05; % in seconds
cfg = ft_definetrial(cfg);

%% artifact rejection
if false %isfile(ARTFILE)
    load(ARTFILE, 'art_eye', 'art_muscle', 'badchannel');
else   
    [art_eye, art_muscle, badchannel] = detect_artifacts(cfg.trl, EEGFILE, ARTFILE);
    save(ARTFILE, 'art_eye', 'art_muscle', 'badchannel');
end
cfg.artfctdef.reject            = 'complete'; 
cfg.artfctdef.eog.artifact      = art_eye; 
cfg.artfctdef.muscle.artifact   = art_muscle;
cfg = ft_rejectartifact(cfg);

%% redefine trials for analysis
cfg = [];
cfg.datafile                = EEGFILE;
cfg.headerfile              = EEGFILE;
cfg.trialfun                = 'ft_trialfun_general'; % default
cfg.trialdef.eventtype      = 'trigger';
cfg.trialdef.eventvalue     = 'TGCI'; 
cfg.trialdef.prestim        = 0.9; % in seconds
cfg.trialdef.poststim       = 0.3; % in seconds
cfg = ft_definetrial(cfg);

%% minimal preprocessing of data
cfg.channel    = 'all';
% re-reference
cfg.reref                   = 'yes';
cfg.refchannel              = 'all'; % a.k.a. average reference 
cfg.implicitref             = 'chan129'; % add reference channel back
% demean
cfg.demean                  = 'yes';
% bandpass filter
cfg.padding                 = 2; 
cfg.bpfilter                = 'yes';
cfg.bpfreq                  = [20 300];
data = ft_preprocessing(cfg);

%% reject artifactual trials
cfg = [];
cfg.artfctdef.reject            = 'complete'; 
cfg.artfctdef.eog.artifact      = art_eye; 
cfg.artfctdef.muscle.artifact   = art_muscle;
data = ft_rejectartifact(cfg, data);

%% rename channels to be consistent with EGI coordinate files
for i = 1:128
    data.label{i} = strrep(data.label{i}, 'chan', 'E');
    data.label{i} = strrep(data.label{i}, 'E00', 'E');
    data.label{i} = strrep(data.label{i}, 'E0', 'E');
end
data.label{129} = 'Cz';


%% remove line noise
% notch filter 60 Hz and harmonics
cfg.dftfilter               = 'yes';
cfg.dftfreq                 = [60, 120, 180, 240];
data = ft_preprocessing(cfg, data);

%% compute scalp current density with surface Laplacian
% while interpolating bad channels
cfg = [];
cfg.elec = ELEC;
cfg.degree = 20; % this is the default for 128 electrodes
cfg.feedback = 'no';
cfg.badchannel     = badchannel;
data = ft_scalpcurrentdensity(cfg, data);


end







