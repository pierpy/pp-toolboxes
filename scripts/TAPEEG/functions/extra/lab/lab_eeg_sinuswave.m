% Calculates eeg signals using sinewave with fixed,random phases
%
% [data,header] = lab_eeg_sinuswave(chans,numtf,fs,freq,randphase,matrix,lag)
%
% chans        = number of channels in result
% numtf        = number of timeframes in result
% fs           = samplingrate
% freq         = frequency for sinewave
% randphase    = true: enabled / false: disabled
% matrix       = matrix with connections (chans x chans)
% lag          = lag for connections in timeframes
%
% written by F. Hatz Vumc 2013

function [data,header] = lab_eeg_sinuswave(chans,numtf,fs,freq,randphase)

if ~exist('randphase','var') | isempty(randphase)
    randphase = false;
end
if ~exist('freq','var') | isempty(freq)
    data = [];
    header = [];
    return
end
if ~exist('numtf','var')
    numtf = 5000;
end
if ~exist('chans','var')
    chans = 1;
end

data = zeros(chans,numtf);
for i = 1:chans    
    data(i,:) = lab_create_sinewave(freq,numtf/fs,fs,randphase);
end

header.numtimeframes = size(data,2);
header.numchannels = size(data,1);
header.numdatachannels = header.numchannels;
header.samplingrate = fs;
header.EEG_file = ['Sinuswave_' num2str(freq) 'Hz_' num2str(size(data,1)) 'x' num2str(size(data,2))];
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

return
