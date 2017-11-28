% Read Cartool .sef
%
% [data,header,cfg] = lab_read_sef(filename,cfg)
%
% Original author: pierre.megevand@gmail.com
% Modified for TAPEEG: F. Hatz 2012

function [data,header,cfg] = lab_read_sef(Filename,cfg,nodata,segment)

if ~exist('segment','var')
    segment = [];
end
if ~exist('nodata','var')
    nodata = false;
end
if ~exist('cfg','var')
    cfg = [];
end

fid=fopen(Filename,'r');

% Read header
header.version=strcat(fread(fid,4,'int8=>char')');
header.numchannels=fread(fid,1,'int32');
header.numauxchannels=fread(fid,1,'int32');
header.numtimeframes=fread(fid,1,'int32');
header.samplingrate=fread(fid,1,'float32');
header.year=fread(fid,1,'int16');
header.month=fread(fid,1,'int16');
header.day=fread(fid,1,'int16');
header.hour=fread(fid,1,'int16');
header.minute=fread(fid,1,'int16');
header.second=fread(fid,1,'int16');
header.millisecond=fread(fid,1,'int16');

% Read channel names
header.channels=fread(fid,[8,header.numchannels],'int8=>char')';

% Read data
if nodata == false
    if ~isempty(segment) & length(segment) > 1
        if segment(1) > header.numtimeframes-1;
            segment(1) = header.numtimeframes-1;
        end
        if segment(2) > header.numtimeframes;
            segment(2) = header.numtimeframes;
        end
        if segment(1) > 1
            fseek(fid,header.numchannels*(segment(1)-1)*4,0);
        end
        data=fread(fid,[header.numchannels,(segment(2) - segment(1) + 1)],'float32');
    else
        data=fread(fid,[header.numchannels,header.numtimeframes],'float32');
    end
else
    data = zeros(header.numchannels,1);
end

% Close file
fclose(fid);
clearvars fid;

% Correct sef channel-names for old script-bug
tmp = unique(header.channels(:,1));
if length(tmp) == 1 & strcmp(tmp,'e')
    header.channels = header.channels(:,2:end);
end
clearvars tmp

return