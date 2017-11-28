% [data,header,cfg] = lab_read_eph(filename,cfg) - read Cartool *.eph
%
% Original author: pierre.megevand@gmail.com
% Modified for TAPEEG: F. Hatz 2012

function [data,header,cfg] = lab_read_eph(filename,cfg)

if ~exist('cfg','var')
    cfg = [];
end

if strcmp(filename(end-3:end),'.eph')==1
    fid=fopen(filename,'rt');
    % read header
    eph=textscan(fid,'%s','delimiter','/n');
    eph=eph{1};
    header.numchannels=sscanf(eph{1},'%f',1);
    header.numtimeframes=sscanf(eph{1},'%*f %f',1);
    header.samplingrate=sscanf(eph{1},'%*f %*f %f',1);
    % prepare for reading data
    formatstring='%f';
    if header.numchannels>1
        for i=1:header.numchannels-1
            formatstring=[formatstring ' %f'];
        end
    end
    
    % read data
    data=zeros(header.numtimeframes,header.numchannels);
    for i=1:header.numtimeframes
        data(i,:)=sscanf(eph{i+1},formatstring);
    end
    fclose(fid);
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
