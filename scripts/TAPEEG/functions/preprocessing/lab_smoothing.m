% Smooth results of signal data with gaussian kernel
%
% [data,header,cfg] = lab_smoothing(data,header,cfg)
%
% Initially distances in loc file are normalized (minimal distance = 1).
% 'Max distance' defines the maximal distance of an electrode from the
% center electrode to be included in the gaussian kernel. 'Weight max
% distance' is the weigth of that electrodes and defines 'sigma' of
% gaussian kernel
%
% Optionally the optimal number of 'Max distance' can be defined by a
% cross-correlation analysis
%
% If called without any input variables, also excel-files are accepted as
% input: The excel file can contain information of multiple variables, the number
% of results values is divided by the number of selected electrode
% locations (locs in electrode file minus excluded channels) to define the
% number of different variables. Remaining variables at the end of the file
% are taken as outcomes and are omitted.
%
% Written by F. Hatz 2011/2012/2014

function [data,header,cfg,skipprocessing] = lab_smoothing(data,header,cfg,nowriting)

disp ('   Smooth data')

skipprocessing = 0;
isexcel = false;
if ~exist('data','var')
    disp('  Select EEG/MEG or Excel-file')
    [data,header,cfg] = lab_read_data;
    if isempty(header) & ~iscell(data)
        disp('   Abort: wrong input data')
        return
    elseif iscell(data)
        isexcel = true;
        if ischar(cfg)
            [header.EEG_file,header.EEG_filepath] = lab_filename(cfg);
            cfg = [];
        end
        disp('  Subjects in row or column')
        Mtranspose = questdlg('Are subjects names in first row or column?','Subjects column/row','Cancel','Column','Row','Row');
        if isempty(Mtranspose) | strcmp(Mtranspose,'Cancel')
            return
        elseif strcmp(Mtranspose,'Column')
            data = data';
        end
        data = lab_correctheader(data);
        [data,cfg] = lab_getstructure(data,cfg);
        if ~isfield(cfg,'clustervars') & cfg.clustervars < 1
            cfg.clustervars = 1;
        end
        cfg.numvars = floor((size(data,1)-1) / cfg.clustervars);
        for i = 1:cfg.clustervars
            tmp = strfind(data{i+1,1},'_');
            if ~isempty(tmp)
                Channels{i,1} = data{i+1,1}(tmp(end)+1:end); %#ok<AGROW>
            else
                Channels{i,1} = data{i+1,1}; %#ok<AGROW>
            end
        end
        data = data(1:cfg.numvars*cfg.clustervars+1,:);
        header.channels = char(Channels);
        header.subjects = data(1,2:end);
        header.vars = data(2:cfg.numvars*cfg.clustervars-1,1);
        clearvars Channels tmp
        try
            datatmp = cell2mat(data(2:end,2:end));
            data = reshape(datatmp,[cfg.clustervars cfg.numvars size(datatmp,2)]);
        catch %#ok<CTCH>
            disp('   Abort: wrong input data')
            return
        end
    end
elseif ~exist('header','var')
    header = lab_create_header(data);
end
if ~exist('cfg','var')
    cfg = [];
end
if ~isfield(cfg,'Output_file') & isfield(header,'EEG_file')
    cfg.Output_file = header.EEG_file;
    cfg.Output_filepath = header.EEG_filepath;
end
if isfield(cfg,'Output_file')
    cfg.Output_fileS = lab_filename(cfg.Output_file);
end
if ~exist('nowriting','var')
    nowriting = false;
end

if ~isfield(cfg,'SMOOTH') | isempty(cfg.SMOOTH)
    [cfg,skipprocessing] = lab_set_smoothing(cfg,header);
    if skipprocessing == 1
        return
    end
end
maxdistance = cfg.SMOOTH.maxdistance;
weightmaxdistance = cfg.SMOOTH.weightmaxdistance;

Idx = 1:size(data,1);
if isexcel == false
    if isfield(header,'numdatachannels')
        Idx = 1:header.numdatachannels;
    end
    if isfield(cfg,'exclude') & ~isempty(cfg.exclude) & (~isfield(cfg.SKIP,'EXCLUDE') | cfg.SKIP.EXCLUDE == false)
        Idx = setdiff(Idx,cfg.exclude);
    end
    [datatmp,headertmp] = lab_reduce_channels(data,header,Idx);
else
    datatmp = data;
    headertmp = header;
end

% Interpolate bad
if isfield(cfg.SMOOTH,'interpolate') & cfg.SMOOTH.interpolate == 1
    disp('   Interpolate bad channels')
    [datatmp,headertmp] = lab_interpolate_bad(datatmp,headertmp);
end

Nchans = size(datatmp,1);
if isfield(headertmp,'locs')
    LOCS = headertmp.locs;
elseif isfield(cfg.SMOOTH,'Locs') & ~isempty(cfg.SMOOTH.Locs)
    LOCS = cfg.SMOOTH.Locs;
    if length(cfg.SMOOTH.Locs.x) ~= size(datatmp,1) & isfield(cfg,'exclude') & ~isempty(cfg.exclude)
        LOCS = lab_reduce_locs(LOCS,cfg.exclude);
    end
    if length(cfg.SMOOTH.Locs.x) ~= size(datatmp,1)
        disp('   Abort: missing/wrong locs-information')
        return
    end
else
    disp('   Abort: missing/wrong locs-information')
    return
end

disp ('  Calculate distance matrix')
distances = lab_distance(LOCS);

disp ('  Normalize distance (minimal distance = 1)')
distances = distances / min(distances(distances > 0));

% create Output_folder
if isfield(cfg,'Output_filepath')
    if ~isfield(cfg.SMOOTH,'folder')
        cfg.SMOOTH.folder = 'Smoothing';
    end
    warning off %#ok<WNOFF>
    mkdir(fullfile(cfg.Output_filepath,cfg.SMOOTH.folder));
    warning on %#ok<WNON>
    Output_filepath = fullfile(cfg.Output_filepath,cfg.SMOOTH.folder);
    Output_fileS = cfg.Output_fileS;
else
    nowriting = true;
end

% Find max distance automatically
if ~isempty(cfg.SMOOTH.AUTO)
    if isfield(cfg,'MAIN') & isfield(cfg.MAIN,'auto') & cfg.MAIN.auto == 1
        maxdistance = lab_test_smoothing(datatmp,LOCS,distances,cfg.SMOOTH.AUTO,fullfile(Output_filepath,Output_fileS),true);
    else
        maxdistance = lab_test_smoothing(datatmp,LOCS,distances,cfg.SMOOTH.AUTO,fullfile(Output_filepath,Output_fileS));
    end
    cfg.SMOOTH.maxdistance = maxdistance;
end

% Calculate neighbours weights
neighbours = zeros(Nchans,Nchans);
sigma =   (maxdistance^2 / (-2 * log(weightmaxdistance / 100)))^0.5;
for i = 1:Nchans;
    for j = 1:Nchans;
        if distances(i,j) <= maxdistance
            neighbours(i,j) = exp(-distances(i,j)^2 / (sigma^2 * 2));
        else
            neighbours(i,j) = 0;
        end
    end
end

% Calculate smoothed data
disp (['  Calculate smoothed data using gaussian kernel (size: ' num2str(maxdistance) ' border-weight: ' num2str(weightmaxdistance) ')'])
for j = 1:size(datatmp,3)
    for i = 1:size(datatmp,1)
        tmp = datatmp(:,i,j)' * neighbours;
        datatmp(:,i,j) = tmp' ./ sum(neighbours,2);
    end
end

% merge data with smoothed data
data(Idx,:,:) = datatmp;

% Write results
if nowriting == false
    disp ('  Write results')
    if isexcel == false
        lab_write_sef(fullfile(Output_filepath,[Output_fileS '_Smooth.sef']),data,header);
    elseif isfield(header,'subjects')
        if isfield(cfg.SMOOTH,'dopower') & cfg.SMOOTH.dopower == true
            dataR = data;
            for i = 1:size(data,3);
                dataR(:,:,i) = dataR(:,:,i) ./ repmat(sum(dataR(:,:,i),2),1,size(dataR,2));
            end
            dataR = reshape(dataR,[size(data,1)*size(data,2) size(data,3)]);
        end
        data = reshape(data,[size(data,1)*size(data,2) size(data,3)]);
        
        xlsout{1,1} = ['C' size(data,1) ' R0'];
        xlsout = [xlsout header.subjects];
        xlsout = cat(1,xlsout,[header.vars num2cell(data)]);
        lab_write_xls(fullfile(Output_filepath,[Output_fileS '_Smooth.xlsx']),xlsout);
        if isfield(cfg.SMOOTH,'dopower') & cfg.SMOOTH.dopower == true
            xlsout{1,1} = ['C' size(data,1) ' R0'];
            xlsout = [xlsout header.subjects];
            xlsout = cat(1,xlsout,[header.vars num2cell(dataR)]);
            lab_write_xls(fullfile(Output_filepath,[Output_fileS '_SmoothR.xlsx']),xlsout);
            
            xlsout{1,1} = ['C' size(data,1) ' R0'];
            xlsout = [xlsout header.subjects];
            xlsout = cat(1,xlsout,[header.vars num2cell(log(data))]);
            lab_write_xls(fullfile(Output_filepath,[Output_fileS '_SmoothLog.xlsx']),xlsout);
            
            xlsout{1,1} = ['C' size(data,1) ' R0'];
            xlsout = [xlsout header.subjects];
            xlsout = cat(1,xlsout,[header.vars num2cell(lab_logit(dataR))]);
            lab_write_xls(fullfile(Output_filepath,[Output_fileS '_SmoothRlogit.xlsx']),xlsout);
        end
    end
end

% Write verbose file (*.vrb)
if exist('Output_filepath','var')
    fid=fopen(fullfile(Output_filepath,[Output_fileS '_smooth.vrb']),'w');
    fprintf(fid,'Smoothing\n');
    fprintf(fid,datestr(now,0));
    fprintf(fid,'\n\n');
    if ~isempty(cfg.SMOOTH.AUTO)
        fprintf(fid,'Find optimal maximal distance: enabled\n');
        fprintf(fid,'\n\n');
    else
        fprintf(fid,'Find optimal maximal distance: disabled\n');
        fprintf(fid,'\n\n');
    end
    fprintf(fid,'Maximal Distance\n');
    fprintf(fid,num2str(maxdistance));
    fprintf(fid,'\n\n');
    fprintf(fid,'Weight Max Distance (percent)\n');
    fprintf(fid,num2str(weightmaxdistance));
    fprintf(fid,'\n\n');
    fprintf(fid,'Sigma\n');
    fprintf(fid,num2str(sigma));
    fprintf(fid,'\n\n');
    fprintf(fid,'Number of Channels\n');
    fprintf(fid,num2str(Nchans));
    fclose(fid);
end

end

