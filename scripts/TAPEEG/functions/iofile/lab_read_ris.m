% Read Cartool .ris
%
% [data,RIS,cfg] = lab_read_ris(filename,cfg)
%
% written by F. Hatz 2012

function [data,RIS,cfg] = lab_read_ris(filename,cfg)

if ~exist('cfg','var')
    cfg = [];
end

if ~exist('filename','var')
    [RIS_file,RIS_filepath]=uigetfile('*.ris','Select is file');
else
    RIS_file=filename;
    tmp=strfind(RIS_file,filesep);
    if ~isempty(tmp)
        RIS_filepath=RIS_file(1:tmp(end));
        RIS_file=RIS_file(tmp(end)+1:end);
    else
        RIS_filepath=[];
    end
    clearvars tmp;
end
fid=fopen(fullfile(RIS_filepath,RIS_file));
if fid > 0
    RIS.version=strcat(fread(fid,4,'int8=>char')');
    if strcmp(RIS.version,'RI01')
        RIS.numsolutionpoints=fread(fid,1,'int32');
        RIS.numtimeframes=fread(fid,1,'int32');
        RIS.samplingrate=fread(fid,1,'float');
        RIS.isscalar=fread(fid,1,'int8');
        if RIS.isscalar == 1
            RIS.data = zeros(RIS.numsolutionpoints,RIS.numtimeframes);
            for i = 1:RIS.numtimeframes
                RIS.data(:,i) = fread(fid,RIS.numsolutionpoints,'float')';
            end
        else
            tmp = zeros(RIS.numsolutionpoints*3,RIS.numtimeframes);
            for i = 1:RIS.numtimeframes
                tmp(:,i) = fread(fid,(RIS.numsolutionpoints * 3),'float')';
            end
            tmp = reshape(tmp,3,RIS.numsolutionpoints,RIS.numtimeframes);
            tmp = permute(tmp,[2 3 1]);
            RIS.x = tmp(:,:,1);
            RIS.y = tmp(:,:,2);
            RIS.z = tmp(:,:,3);
            RIS.data = (RIS.x.^2 + RIS.y.^2 + RIS.z.^2).^0.5;
        end
    end
    fclose(fid);
    RIS.numchannels = RIS.numsolutionpoints;
    RIS.channels = num2str((1:RIS.numchannels)');
    data = RIS.data;
else
    RIS = [];
end
return