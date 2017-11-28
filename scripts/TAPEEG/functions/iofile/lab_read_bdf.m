% Reads .bdf NXE continous datafile

function [data,header,cfg] = lab_read_bdf(Filename,cfg)

if exist(Filename,'file')
    
hdr = read_biosemi_bdf(Filename);
header.samplingrate = hdr.Fs;
header.numchannels = hdr.nChans;
header.numtimeframes = hdr.nSamples * hdr.nTrials;
header.channels = char(hdr.labels(:));
data = read_biosemi_bdf(Filename,hdr,1,header.numtimeframes,1:header.numchannels);

end