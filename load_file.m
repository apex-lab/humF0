function TFR = load_file(i, f0) 

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

%% define trials
cfg = [];
cfg.datafile                = EEGFILE;
cfg.headerfile              = EEGFILE;
cfg.trialfun                = 'ft_trialfun_general'; % default
cfg.trialdef.eventtype      = 'trigger';
cfg.trialdef.eventvalue     = 'TGCI'; 
cfg.trialdef.prestim        = 0.7; % in seconds
cfg.trialdef.poststim       = 0.0; % in seconds
cfg = ft_definetrial(cfg);

%% artifact rejection
if isfile(ARTFILE)
    load(ARTFILE, 'art_eye', 'art_muscle', 'badchannel');
else   
    [art_eye, art_muscle, badchannel] = detect_artifacts(cfg.trl, EEGFILE);
    save(ARTFILE, 'art_eye', 'art_muscle', 'badchannel');
end
cfg.artfctdef.reject            = 'complete'; 
cfg.artfctdef.eog.artifact      = art_eye; 
cfg.artfctdef.muscle.artifact   = art_muscle;
cfg = ft_rejectartifact(cfg);

%% minimal preprocessing of artifact free data
cfg.channel    = 'all';
% re-reference
cfg.reref                   = 'yes';
cfg.refchannel              = 'all'; % a.k.a. average reference 
cfg.implicitref             = 'chan129'; % add reference channel back
% bandpass filter
cfg.padding                 = 2; 
cfg.bpfilter                = 'yes';
cfg.bpfreq                  = [20 300];
data = ft_preprocessing(cfg);

%% repair bad channels
cfg         = [];
cfg.elec    = ELEC;
cfg.method  = 'distance';
neighbours  = ft_prepare_neighbours(cfg, data);

cfg = [];
cfg.badchannel     = badchannel;
cfg.method         = 'weighted';
cfg.neighbours     = neighbours;
data = ft_channelrepair(cfg, data);

%% now detrend the data
cfg             = [];
cfg.channel     = 'all';
cfg.demean      = 'yes';
cfg.polyremoval = 'yes';
cfg.polyorder   = 1; 
data = ft_preprocessing(cfg, data);

%% robustly re-reference detrended data (de Cheveigne & Arzounian, 2018) 
cfg                         = [];
cfg.channel                 = 'all';
% re-reference since average is no longer the same after interpolation
cfg.reref                   = 'yes';
cfg.refchannel              = 'all'; 
data = ft_preprocessing(cfg, data);

%% remove line noise
% notch filter 60 Hz and harmonics
cfg.dftfilter               = 'yes';
cfg.dftfreq                 = [60, 120, 180, 240];
data = ft_preprocessing(cfg, data);


%% time-frequency analysis (wavelet)
cfg                  = [];
cfg.channel          = 'all';
cfg.method           = 'wavelet';
cfg.width            = round(0.06*f0*pi); % wavelet length of ~60ms for F0
cfg.output           = 'pow';
f0 = round(f0);
cfg.foi              = (f0 - 60):1:(f0 + 60);
cfg.toi              = -0.6:0.01:0;   
cfg.keeptrials       = 'no';
TFR = ft_freqanalysis(cfg, data);  
TFR.freq = -60:1:60; % normalize by F0 for later averaging across S's


end







