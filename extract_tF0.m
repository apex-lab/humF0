function tF0 = extract_tF0(subNum)

    Fs = 48000;

    %%  load file
    s = num2str(subNum);
    while length(s) < 3
        s = strcat('0', s);
    end
    try
        fname = strcat("/home/john/Documents/EGG/sub-P", s, "/ses-S001/eeg/sub-P", s, "_ses-S001_task-slowHum_run-001_eeg.mat");
        load(fname, 'xdf');
    catch
        fname = strcat("/home/john/Documents/EGG/sub-P", s, "/ses-S001/eeg/sub-P", s, "_ses-S001_task-slowHum_run-002_eeg.mat");
        load(fname, 'xdf');
    end
     
        
    %% loop through trials and extract f0 from 100-200 ms post GCI
    % extract time of glottis closure instants 
    for i = 1:length(xdf)
        if (xdf{i}.info.name == "glottis_closure_instants")
            stream = i;
        end
    end
    gci = xdf{stream}.time_stamps;
    % extract EGG trace
    for i = 1:length(xdf)
        if (xdf{i}.info.name == "AudioCaptureWin")
            stream = i;
        end
    end
    egg = xdf{stream}.time_series(1,:);
    t = xdf{stream}.time_stamps;
    % get indices for glottis closures in EGG trace
    m = 1:length(gci);
    for i = 1:length(gci)
        m(i) = find(t > gci(i), 1);
    end
    sum_f0 = zeros(35,1); % running sum
    for i = 1:length(m)
        % get 00 ms to 500 ms of each trial
        trial = egg(m(i) + 0.0*Fs:m(i) + 0.5*Fs); 
        
        % and do an FFT
%         L = length(trial);    % Length of signal
%         Y = fft(trial);
%         P2 = abs(Y/L);
%         P1 = P2(1:L/2+1);   
        %P1(2:end-1) = 2*P1(2:end-1);
        %f = Fs*(0:(L/2))/L; % frequency in Hz
        
        % peak power of EGG is fundamental frequency
        %[~, idx] = max(P1);
        %f0 = f(idx);
        
        [p, loc] = pitch(trial', Fs);
        %f0 = mean(p);
        if (i == 1)
            sum_f0 = p; 
        else
            sum_f0 = sum_f0 + p; 
        end
    end

    %% average terminal F0 and then return
    F0 = sum_f0/length(m);
    plot(loc/Fs, F0)
    [~, tF0, ~] = ginput(1); % click the terminal F0 on figure
end