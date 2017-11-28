% Calculates distance-matrices of localisation-information in header.locs
%
% [result,cfg] = lab_calculate_distance(data,header,cfg)
% [result,cfg] = lab_calculate_distance(header)
%
% written by F. Hatz 2013

function [result,cfg] = lab_calculate_distance(data,header,cfg)
disp('Calculate distance matrix')

if isstruct(data) & isfield(data,'locs')
    header = data;
end
if ~exist('cfg','var')
    cfg = [];
end
if isfield(cfg,'SKIP') & isfield(cfg.SKIP,'DISTANCE') & cfg.SKIP.DISTANCE == true;
    return
end

if ~isfield(header,'EEG_file')
   if isfield(cfg,'Output_file')
       header.EEG_file = cfg.Output_file;
   else
       header.EEG_file = '';
   end
end
if ~isfield(cfg,'Output_file') & isfield(header,'EEG_filepath')
    cfg.Output_filepath = header.EEG_filepath;
    cfg.Output_file = header.EEG_file;
end
if isfield(cfg,'Output_file')
    [~,~,~,cfg.Output_fileS] = lab_filename(cfg.Output_file);
end

result.distance = lab_distance(header.locs);
if max(result.distance(:)) < 1
    result.distance = result.distance * 1000;
end

% save Result
if isfield(cfg,'Output_filepath') & ~isempty(cfg.Output_filepath) & ~exist('dosingle','var')
    warning off %#ok<WNOFF>
    mkdir(fullfile(cfg.Output_filepath,'Distance'));
    Distance_filepath = fullfile(cfg.Output_filepath,'Distance');
    if isfield(cfg,'patient')
        Distance_file = [cfg.patient '_Distance.mat'];
    else
        Distance_file = 'Distance-mat';
        cfg.patient = [];
    end
    if ~isfield(cfg,'listold')
        if exist(fullfile(Distance_filepath,Distance_file),'file')
            disp('     delete Distance.mat from previous run')
            delete(fullfile(Distance_filepath,Distance_file));
        end
        cfg.listold = cellstr(fullfile(Distance_filepath,Distance_file));
    elseif ~strcmp(cfg.listold,fullfile(Distance_filepath,Distance_file))
        if exist(fullfile(Distance_filepath,Distance_file),'file')
            disp('     delete Distance.mat from previous run')
            delete(fullfile(Distance_filepath,Distance_file));
        end
        cfg.listold = [cfg.listold cellstr(fullfile(Distance_filepath,Distance_file))];
    end
    warning on %#ok<WNON>
    
    if exist(fullfile(Distance_filepath,Distance_file),'file')
        result2 = result;
        clearvars result
        try %#ok<TRYNC>
            load(fullfile(Distance_filepath,Distance_file));
        end
        if ~exist('result','var') | ~isfield(result,'file')
            result.file = [];
        end
        if ~isfield(result,'distance')
            result.distance = [];
        end
        if isempty(header.EEG_file)
            result.file = [result.file cellstr(num2str(size(result.distance,3)))];
        else
            result.file = [result.file cellstr(header.EEG_file)];
        end
        result.distance = cat(3,result.distance,result2.distance);
        clearvars result2
    elseif isempty(header.EEG_file)
        result.file = cellstr(num2str(size(result.distance,3)));
    else
        result.file = cellstr(header.EEG_file);
    end
    patient = cfg.patient; %#ok<NASGU>
    for nmatrix = 1:size(result.distance,3)
        distancefile = fullfile(Distance_filepath,[Distance_file(1:end-4) '_' num2str(nmatrix,'%03d') '.txt']);
        dlmwrite(distancefile,result.distance(:,:,nmatrix),'delimiter','\t','precision', 6);
    end
    save(fullfile(Distance_filepath,Distance_file),'result','patient');
end