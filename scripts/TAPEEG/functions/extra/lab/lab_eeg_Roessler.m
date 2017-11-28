function [data,header] = lab_eeg_Roessler(nchans,length,fs,settings)

if ~exist('settings','var')
    settings = [];
end
if ~isfield(settings,'a') | isempty(settings.a)
    settings.a = 0.2;
end
if ~isfield(settings,'b') | isempty(settings.b)
    settings.b = 0.4;
end
if ~isfield(settings,'c') | isempty(settings.c)
    settings.c = 5.7;
end
if ~isfield(settings,'noise') | isempty(settings.noise)
    settings.noise = 0;
end

data = zeros(nchans,length*fs);
for i = 1:nchans
    [x,y,z] = rossler(length*fs + 500,settings.noise,settings.a,settings.b,settings.c);
    signal = (x.^2 + y.^2 + (z/4).^2).^0.5;
    signal = signal(501:end);
    signal = signal - mean(signal);
    data(i,:) = signal(:)';
end

header.numtimeframes = size(data,2);
header.numchannels = size(data,1);
header.numdatachannels = header.numchannels;
header.samplingrate = fs;
header.EEG_file = ['Roessler-model' num2str(size(data,1)) 'x' num2str(size(data,2))];
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