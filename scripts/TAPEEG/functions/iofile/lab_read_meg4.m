% [data,header,cfg] = lab_read_meg4(filename,cfg) - read CTF-MEG
%
% written by F. Hatz 2012

function [data,header,cfg] = lab_read_meg4(filename,cfg)

if ~exist('cfg','var')
    cfg = [];
end
if ~exist('filename','var')
    [filename,filepath] = uigetfile('*.fif;*.fiff','Select file');
    filename = fullfile(filepath,filename);
end

warning off %#ok<WNOFF>
hdr = ft_read_header(filename);
data = ft_read_data(filename);
warning on %#ok<WNON>
if size(data,3) > 1
    for i = 1:size(data,3)
        events.POS(1,i) = int64((i-1)*size(data,2)+1);
        events.DUR(1,i) = int64(1);
        events.OFF(1,i) = int64(0);
        events.TYP{1,i} = ['Trial' num2str(i,'%03d')];
    end
    data = reshape(data,size(data,1),size(data,2)*size(data,3));
else
    events = [];
end
if strcmp(hdr.label{1,1},'STIM')
    eventstmp = data(1,:);
    tmp1 = unique(eventstmp);
    events2.POS = [];
    events2.DUR = [];
    events2.OFF = [];
    events2.TYP = {};
    for i = 1:length(tmp1)
        tmp2 = find(eventstmp == tmp1(i));
        events2.POS = [events2.POS int64(tmp2)];
        events2.DUR = [events2.DUR int64(ones(1,length(tmp2)))];
        events2.OFF = [events2.OFF int64(zeros(1,length(tmp2)))];
        events2.TYP = [events2.TYP repmat(cellstr(['Trig' num2str(i,'%03d')]),[1 length(tmp2)])];
    end
    events = lab_mix_markers(events,events2);
    data = data(2:end,:);
    hdr.label = hdr.label(2:end,:);
    hdr.nChans = hdr.nChans - 1;
end

header.year=0;
header.month=0;
header.day=0;
header.hour=0;
header.minute=0;
header.second=0;
header.millisecond=0;
header.numauxchannels = 0;
header.numchannels = size(data,1);
header.numtimeframes = size(data,2);
header.channels = char(hdr.label);
header.samplingrate = hdr.Fs;

header.orig = hdr;

% create locs
grad = ft_convert_units(hdr.grad,'mm');
header.locs.x = grad.chanpos(:,1)';
header.locs.y = grad.chanpos(:,2)';
header.locs.z = grad.chanpos(:,3)';
header.locs.labels = grad.label';
header.locs.grad = grad;
header.locs = lab_locs2sph(header.locs);

header.datatype = 'meg';
header.ref_chan = 'none';

% Define data channels
header.numauxchannels = header.numchannels - size(header.locs.x,2);
header.numdatachannels = header.numchannels - header.numauxchannels;

% Set events
header.events = events;

return






