% Read Netstation raw export
%
% [data,header,cfg] = lab_read_raw(filename,cfg)
%
% written by F. Hatz 2012

function [data,header,cfg] = lab_read_raw(filename,cfg,nodata,segment)

if ~exist('segment','var')
    segment = [];
end
if ~exist('nodata','var')
    nodata = false;
end
if ~exist('cfg','var')
    cfg = [];
end

if nodata == false
    [headertmp,~] = lab_readegi(filename,0);
    if ~isempty(segment) & length(segment) > 1
        segmented = 0;
        switch(headertmp.version),
            case {2,4,6}
                segmented = 0;
            case {3,5,7}
                segmented = 1;
        end
        if segmented == 0
            [headertmp,data,EventData] = lab_readegi(filename,segment(1):segment(2));
        else
            [headertmp,data,EventData] = lab_readegi(filename);
            if segment(2) > size(data,2)
                segment(2) = size(data,2);
            end
            data = data(:,segment(1):segment(2));
            EventData = EventData(:,segment(1):segment(2));
        end
    else
        [headertmp,data,EventData] = lab_readegi(filename);
    end
else
    [headertmp,data,EventData] = lab_readegi(filename,0);
end
header.samplingrate = headertmp.samp_rate;
header.numchannels = headertmp.nchan;
data = data(1:header.numchannels,:);
header.version=headertmp.version;
header.numauxchannels=0;
if nodata == true
    header.numtimeframes = headertmp.samples;
else
    header.numtimeframes = size(data,2);
end
header.year=0;
header.month=0;
header.day=0;
header.hour=0;
header.minute=0;
header.second=0;
header.millisecond=0;
events.POS = [];
events.DUR = [];
events.OFF = [];
events.TYP = [];
if find(EventData > 0) > 0
    tmp = zeros(1,size(EventData,2));
    for i = 1:size(headertmp.eventcode,1)
        tmp(1,EventData(i,:) == 1) = i;
    end
    tmp2 = find(tmp > 0);
    tmp2(2,:) = tmp(1,tmp2);
    for i = 1:size(tmp2,2)
        events.POS = [events.POS int64(tmp2(1,i))];
        events.DUR = [events.DUR int64(1)];
        events.OFF = [events.OFF int64(0)];
        events.TYP = [events.TYP cellstr(headertmp.eventcode(tmp2(2,i),:))];
    end
    header.events = events;
end
clearvars headertmp CatIndex;
if header.numchannels == 256
    data(end+1,:) = zeros(1,size(data,2));
    header.numchannels = 257;
end
header.channels = num2str((1:header.numchannels)');
