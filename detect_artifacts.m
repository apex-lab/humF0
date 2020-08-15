function [art_eye, art_muscle, badchannel] = detect_artifacts(trl, datafile, artfile)
% identifies artifacts by (1) reading the trials with padding from disk,
% (2) filtering the data to the frequency range target artifacts are
% known to occur in, (3) z-transforming the Hilbert envelope of the data,
% and (4) rejecting when the z-value passes a generous threshold, since
% artifacts can easily pass such a threshold. In practice, you will want to
% enable the lines of code that make the z-threshold interactive so you
% can make sure that all artifacts are being rejected and not too much
% else, but it is nice to automate the end procedure for reproducibility.

    %% define EOG channels
    EOG = [126 8 127 25 17 1 32]; % Jia and Tyler (2019)

    %% manually remove bad channels
    % use "identify" button to get names of channels you see have problems
    if ~isfile(artfile)
    cfg = [];
    cfg.dataset = datafile;
    cfg.datafile = datafile;
    cfg.headerfile = datafile; 
    cfg.continuous = 'yes';
    cfg.channel       = 'all';
    cfg.viewmode      = 'vertical';
    cfg.blocksize     = 10; % time window to browse
    cfg.preproc.dftfilter = 'yes';
    cfg.preproc.dftfreq = [60, 120, 180, 240];
    cfg.preproc.demean = 'yes';
    cfg.preproc.reref = 'yes';
    cfg.preproc.refchannel = 'all';
    cfg.ylim = [-15 15];
    ft_databrowser(cfg);
        
    
    % and type in those names below 
    badchannel  = input('write badchannels in matlab list syntax: ');
    
    else % put them back in the format we use in this function
        load(artfile, 'badchannel');
        new = 1:length(badchannel);
        for i = 1:length(badchannel)
            num = erase(badchannel{i},'E');
            new(i) = str2num(num);
        end
        badchannel = new;
    end
    
    % we'll exclude them from artifact rejection below and later
    % we'll interpolate them in the main analysis
    idx = 1:129;
    good_idx = idx(~ismember(idx, badchannel));
    
    % we also want only the good eog channels
    good_eog = EOG(~ismember(EOG, badchannel));
    
    %% load continuous data (memory intensive)
    % the only purpose of this is to rereference in case the reference 
    % channel has noise in it, otherwise it would be better not
    % to load the full dataset into memory (in which case you just remove
    % the data argument from ft_artifact_zvalue() below and uncomment the
    % datafile lines of cfg)
    cfg = [];
    cfg.dataset     = datafile;
    cfg.preproc.reref = 'yes';
    cfg.preproc.refchannel = 'all';
    data        = ft_preprocessing(cfg);
    
    
    %% identify eye artifacts
    cfg = [];
    cfg.trl = trl; % trial definition
    %cfg.datafile = datafile;
    %cfg.headerfile = datafile; 
    cfg.artfct.reject = 'complete'; % reject full trial of artifact
    cfg.continuous = 'yes';
    % channel selection and padding
    cfg.artfctdef.zvalue.channel = good_eog;
    cfg.artfctdef.zvalue.cutoff = 10;
    cfg.artfctdef.zvalue.trlpadding = 0;
    cfg.artfctdef.zvalue.artpadding = 0;
    cfg.artfctdef.zvalue.fltpadding = 1;
    % now set filtering parameters to emphasize artifact of interest (eye)
    cfg.artfctdef.zvalue.bpfilter = 'yes';
    %cfg.artfctdef.zvalue.bpfilttype = 'but'; % butterworth filter
    cfg.artfctdef.zvalue.bpfiltord = 2;
    cfg.artfctdef.zvalue.bpfreq = [2 15];
    % set to use signal envelope instead of raw signal
    cfg.artfctdef.zvalue.hilbert = 'yes'; 
    %cfg.artfctdef.zvalue.interactive='yes';
    [~, art_eye] = ft_artifact_zvalue(cfg, data);
    
    %% identify muscle artifacts
    cfg = [];
    cfg.trl = trl; 
    %cfg.datafile = datafile;
    %cfg.headerfile = datafile; 
    cfg.artfct.reject = 'complete'; 
    cfg.continuous = 'yes';
    % channel selection and padding
    cfg.artfctdef.zvalue.channel = good_idx; 
    cfg.artfctdef.zvalue.cutoff = 10;
    cfg.artfctdef.zvalue.trlpadding = 0;
    cfg.artfctdef.zvalue.artpadding = 0;
    cfg.artfctdef.zvalue.fltpadding = 1;
    % filter parameters
    cfg.artfctdef.zvalue.bpfilter = 'yes';
    cfg.artfctdef.zvalue.bpfilttype = 'but'; % butterworth filter
    cfg.artfctdef.zvalue.bpfreq = [110 140];
    cfg.artfctdef.zvalue.dftfilter = 'yes';
    cfg.artfctdef.zvalue.dftfreq = [60, 120, 180, 240];
    % set to use signal envelope instead of raw signal
    cfg.artfctdef.zvalue.hilbert = 'yes'; 
    cfg.artfctdef.zvalue.boxcar = 0.2;
    %cfg.artfctdef.zvalue.interactive='yes';
    [~, art_muscle] = ft_artifact_zvalue(cfg, data);
    
    %% put bad channels in fieldtrip readable format
    k = length(badchannel);
    b = cell(1, k);
    for i = 1:k
        b{i} = strcat('E', num2str(badchannel(i)));
    end
    badchannel = b;

end