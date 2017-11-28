function [data,header] = lab_create_eeg_phaseresetting(chans,numtf,fs,settings)

if ~exist('fs','var') | isempty(fs)
    fs = 250;
end
if ~exist('numtf','var')
    numtf = 4000;
end

%general parameters of the signal
if isfield(settings,'sweepdur')
    FRAMES = settings.sweepdur;
else
    FRAMES = 200;
end
TRIALS = ceil((numtf)/200);
SRATE = fs;
%parameters of NE
if isfield(settings,'nfreq')
    NEFREQ = settings.nfreq;
    NEAMP = settings.namp;
    NEPOS = settings.npos;
else
    NEAMP = -24;
    NEFREQ = 5;
    NEPOS = 115;
end
%parameters of PE
if isfield(settings,'nfreq')
    PEFREQ = settings.pfreq;
    PEAMP = settings.pamp;
    PEPOS = settings.ppos;
else
    PEAMP = 11;
    PEFREQ = 1;
    PEPOS = 150;
end
%temporal jitter of peaks and amplitude of noise
if isfield(settings,'tjitter')
    TJITTER = settings.tjitter;
else
    TJITTER = 8;
end
if isfield(settings,'noiseamp')
    NOISEAMP = settings.noiseamp;
else
    NOISEAMP = 10;
end

disp ('Please wait while generating the data');

load dipole

data = [];
for i = 1:ceil(chans/31)
    ne = NEAMP * peak (FRAMES, TRIALS, SRATE, NEFREQ, NEPOS, TJITTER);
    pe = PEAMP * peak (FRAMES, TRIALS, SRATE, PEFREQ, PEPOS, TJITTER, 1);
    datatmp = dipole (:,1) * ne + dipole (:,4) * pe; %#ok<NODEF>
    data = cat(1,data,datatmp);
end
data = data(1:chans,:);

for ch = 1:size(data,1)
 disp(sprintf('Generating noise for channel %d', ch)); %#ok<DSPS>
 data(ch,:) = data(ch,:) + NOISEAMP * noise (FRAMES, TRIALS, SRATE);
end

header.numtimeframes = size(data,2);
header.numchannels = size(data,1);
header.numdatachannels = header.numchannels;
header.samplingrate = fs;
header.EEG_file = ['PhaseResetting_' num2str(size(data,1)) 'x' num2str(size(data,2))];
header.EEG_filepath = [];
header.channels = num2str((1:size(data,1))');
header.numauxchannels = 0;
tmp = clock;
header.year = tmp(1);
header.month = tmp(2);
header.day = tmp(3);
header.hour = tmp(4);
header.minute = tmp(5);
header.second = floor(tmp(6));
header.millisecond = (tmp(6) - floor(tmp(6)))*1000;
header.ref_chan = 'none';