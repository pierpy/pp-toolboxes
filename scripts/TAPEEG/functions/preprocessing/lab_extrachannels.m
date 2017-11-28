% Set extra channels for eeg/meg-files
%
% [data,header,cfg,skipprocessing] = lab_extrachannels(data,header,cfg)
%
% Written by F. Hatz 2013

function [data,header,cfg,skipprocessing] = lab_extrachannels(data,header,cfg)

skipprocessing = 0;

if ~exist('header','var')
    header =[];
end
if ~exist('data','var')
    data =[];
end
if ~exist('cfg','var')
    cfg = [];
end
if ~isfield(header,'numchannels')
    header.numchannels = size(data,1);
end
if ~isempty(data) & size(data,1) ~= header.numchannels
    header.numchannels = size(data,1);
end
if ~isfield(header,'numdatachannels')
    if isfield(header,'locs') & ~isempty(header.locs) & length(header.locs.x) < header.numchannels
        header.numdatachannels = length(header.locs.x);
    elseif isfield(header,'numchannels')
        header.numdatachannels = header.numchannels;
    end
end

if ~isfield(cfg,'EXTRA') | ~isfield(cfg.EXTRA,'numdatachans') | ~isfield(cfg.EXTRA,'ref_chan') | ...
       ~isfield(cfg.EXTRA,'auxmethod') | ~isfield(cfg.EXTRA,'ecg_ch') | ~isfield(cfg.EXTRA,'computegfp')
    [cfg,skipprocessing] = lab_set_extrachannels(cfg,header,data);
    if skipprocessing == 1
        return
    end
end

% apply aux channel info
if ~isempty(data) & isstruct(header) & isfield(cfg.EXTRA,'auxmethod') & isfield(cfg.EXTRA,'AUX')
    if strcmp(cfg.EXTRA.auxmethod,'channels')
        if isempty(cfg.EXTRA.AUX) | (isfield(cfg.EXTRA.AUX,'channels') & max(cfg.EXTRA.AUX.channels) > size(data,1))
            cfg.EXTRA.auxmethod = 'automatic';
            cfg.EXTRA.AUX = [];
            disp('AUX channels not in data-range, AUX method set to automatic')
        end
    end
    if strcmp(cfg.EXTRA.auxmethod,'none')
        if isfield(header,'numdatachannels')
            [data,header] = lab_reduce_channels(data,header,(1:header.numdatachannels));
        else
            header.numdatachannels = header.numchannels;
        end
        header.numauxchannels = 0;
    elseif strcmp(cfg.EXTRA.auxmethod,'channels')
        [auxintersect,~,auxIdx] =  intersect((1:header.numdatachannels),cfg.EXTRA.AUX.channels);
        if max(setdiff((1:header.numdatachannels),auxintersect)) < min(auxintersect)
            header.numauxchannels = header.numauxchannels + length(intersect((1:header.numdatachannels),cfg.EXTRA.AUX.channels));
            header.numdatachannels = header.numchannels - header.numauxchannels;
        elseif ~isempty(auxintersect)
            auxchannels = data(auxintersect,:);
            data = cat(1,data,auxchannels);
            header.numauxchannels = header.numauxchannels + size(auxchannels,1);
            header.numchannels = size(data,1);
            clearvars includechannels auxchannels
            tmp = cellstr(header.channels);
            if isfield(cfg.EXTRA.AUX,'labels') & ~isempty(cfg.EXTRA.AUX.labels)
                for i = auxIdx
                    tmp = cat(1,tmp,cfg.EXTRA.AUX.labels(i));
                end
            else
                for i = 1:length(auxintersect)
                    tmp = cat(1,tmp,cellstr(['AUX' num2str(i)]));
                end
            end
            header.channels = char(tmp);
            clearvars tmp i
        end
    end
end

% control conflict number data channels
if isfield(header,'numchannels') & cfg.EXTRA.numdatachans < header.numchannels;
    header.numdatachannels = cfg.EXTRA.numdatachans;
    header.numauxchannels = header.numchannels - cfg.EXTRA.numdatachans;
end

% control conflict and set ref_chan
if isfield(header,'ref_chan')
    if ~isempty(header.ref_chan)
        nonequal = 1;
        if ischar(cfg.EXTRA.ref_chan) & ischar(header.ref_chan) & ...
                strcmp(cfg.EXTRA.ref_chan,header.ref_chan)
            nonequal = 0;
        elseif isnumeric(cfg.EXTRA.ref_chan) & isnumeric(header.ref_chan) & ...
                length(cfg.EXTRA.ref_chan) == length(header.ref_chan) & ...
                min(cfg.EXTRA.ref_chan == header.ref_chan)
            nonequal = 0;
        end
        if nonequal == 1
            if ~isfield(cfg.EXTRA,'REFoverwrite')
                if isfield(cfg,'MAIN') & isfield(cfg.MAIN,'auto') & cfg.MAIN.auto == 1
                    cfg.EXTRA.REFoverwrite = true;
                else
                    if ~isfield(cfg,'EXTRA') | ~isfield(cfg.EXTRA,'REFoverwrite')
                        defanswer = 'Yes';
                    elseif cfg.EXTRA.REFoverwrite == true;
                        defanswer = 'Yes';
                    else
                        defanswer = 'No';
                    end
                    answer = questdlg('Overwrite file-reference channel?','Overwrite ref-chan','No','Yes',defanswer);
                    pause(0.1)
                    if strcmp(answer,'Yes')
                        cfg.EXTRA.REFoverwrite = true;
                    else
                        cfg.EXTRA.REFoverwrite = false;
                    end
                    clearvars answer defanswer
                end
            end
            if cfg.EXTRA.REFoverwrite == false
                cfg.EXTRA.ref_chan = header.ref_chan;
            end
        end
    end
    if isnumeric(cfg.EXTRA.ref_chan) & isfield(header,'numchannels') & ...
            header.numchannels < max(cfg.EXTRA.ref_chan)
        cfg.EXTRA = rmfield(cfg.EXTRA,'ref_chan');
        disp('Abort: Reference channel not in data-range')
        skipprocessing = 1;
        return
    elseif min(cfg.EXTRA.ref_chan) <= 0
        header.ref_chan = 'none';
    elseif ~isfield(cfg,'ADDREF') | ~isfield(cfg.ADDREF,'name') | isempty(cfg.ADDREF.name)
        header.ref_chan = cfg.EXTRA.ref_chan;
    end
end

% control and correct for data in ref_chan (if numeric)
if isfield(header,'ref_chan') & isnumeric(header.ref_chan) & ...
        ~isempty(header.ref_chan) & ~isempty(data)
    dataref = zeros(length(header.ref_chan),size(data,2));
    for i = 1:length(header.ref_chan)
        dataref(i,:) = data(header.ref_chan(i),:);
    end
    if max(abs(dataref(:))) > 0
        if ~isfield(cfg.EXTRA,'REFcorrect')
            if isfield(cfg,'MAIN') & isfield(cfg.MAIN,'auto') & cfg.MAIN.auto == 1
                cfg.EXTRA.REFcorrect = false;
            else
                answer = questdlg('Reference channel not empty, substract reference channel from other channels?','Reset ref-chan','No','Yes','No');
                pause(0.1)
                if strcmp(answer,'Yes')
                    cfg.EXTRA.REFcorrect = true;
                else
                    cfg.EXTRA.REFcorrect = false;
                end
            end
        end
        if cfg.EXTRA.REFcorrect == false
            disp('    Warning: Reference channel(s) contains singal, reference channel is set to zero')
        else
            disp('    Warning: Reference channel(s) contains singal, substract reference channel from other channels')
            if length(header.ref_chan) > 1
                reference = mean(data(header.ref_chan,:),1);
            else
                reference = data(header.ref_chan,:);
            end
            for i = 1:header.numdatachannels
                data(i,:) = data(i,:) - reference;
            end
        end
        data(header.ref_chan,:) = 0;
    end
end

if isfield(cfg.EXTRA,'computegfp') & cfg.EXTRA.computegfp == true
    [data,header] = lab_compute_gfp(data,header);
end

if (~isfield(header,'ecg_ch') | isempty(header.ecg_ch)) & isfield(cfg.EXTRA,'ecg_ch') & ...
        cfg.EXTRA.ecg_ch > 0 & cfg.EXTRA.ecg_ch <= size(data,1)
    header.ecg_ch = cfg.EXTRA.ecg_ch;
end

return