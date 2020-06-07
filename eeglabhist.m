% converts .mat files containing LSL streams into EEGLab .set files
for i = 7:21

EEG.etc.eeglabvers = 'development head'; % this tracks which version of EEGLAB is being used, you may ignore it
s = num2str(i);
while length(s) < 3
    s = strcat('0', s);
end
try
    fname = strcat("/home/john/Documents/EGG/sub-P", s, "/ses-S001/eeg/sub-P", s, "_ses-S001_task-slowHum_run-001_eeg.mat");
    EEG = pop_loadxdf(fname, 'streamtype', 'EEG', 'exclude_markerstreams', {});
catch
    fname = strcat("/home/john/Documents/EGG/sub-P", s, "/ses-S001/eeg/sub-P", s, "_ses-S001_task-slowHum_run-002_eeg.mat");
    EEG = pop_loadxdf(fname, 'streamtype', 'EEG', 'exclude_markerstreams', {});
end
EEG.setname= strcat('P', s);
EEG = eeg_checkset( EEG );
EEG = pop_saveset( EEG, 'filename',strcat('P', s, '.set'),'filepath','/home/john/Documents/EGG/slowHum/');
EEG = eeg_checkset( EEG );

end
