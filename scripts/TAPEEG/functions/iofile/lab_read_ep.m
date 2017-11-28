% [data,header,cfg] = lab_read_ep(filename,cfg) - read Cartool *.ep
% written by F. Hatz 2012

function [data,header,cfg] = lab_read_ep(filename,cfg)

if ~exist('cfg','var')
    cfg = [];
end

if strcmp(filename(end-2:end),'.ep')==1
    fid=fopen(filename,'rt');
    eph=textscan(fid,'%s','delimiter','/n');
    eph=eph{1};
    for i = 1:size(eph,1)
        data(i,:)=sscanf(eph{i},' %f');
    end
    fclose(fid);
    
    header.numtimeframes=size(data,1);
    header.numchannels=size(data,2);
    header.samplingrate=[];
    data = data';
    
    % set missing header information
    header.channels=char(32*ones(header.numchannels, length(header.numchannels)+1));
    for i = 1:header.numchannels
        header.channels(i,1:(length(num2str(i))+1)) = ['E' num2str(i)];
    end
    clearvars i
    header.year=0;
    header.month=0;
    header.day=0;
    header.hour=0;
    header.minute=0;
    header.second=0;
    header.millisecond=0;
    header.numauxchannels = 0;
else
    error('incorrect file type');
end


