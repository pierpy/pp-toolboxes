% Reads .nxe NXE continous datafile

function [data,header,cfg] = lab_read_nxe(Filename,cfg)

if exist(Filename,'file')
    
hdr = read_nexstim_nxe(Filename);
header.samplingrate = hdr.Fs;
header.numchannels = hdr.nChans;
header.numtimeframes = hdr.nSamples * hdr.nTrials;
header.channels = char(hdr.labels(:));
data = read_nexstim_nxe(Filename,1,header.numtimeframes,1:header.numchannels);

end