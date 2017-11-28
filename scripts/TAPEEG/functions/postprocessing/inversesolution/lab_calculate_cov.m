% Calculate covariance matrix for inverse solutions
%
% [header,settings] = lab_calculate_cov(data,header,settings,cfg)
%
% Optionally covariance can be calculated based on all files found in same
% folder as processed file, marker-episodes can be excluded
%
% written by F. Hatz 2013

function [header,settings] = lab_calculate_cov(data,header,settings,cfg)

if isfield(header,'cov')
    return
end

disp('    Calculate covariance matrix')
if ~exist('settings','var') | ~isfield(settings,'COV') | ~isfield(settings.COV,'covmethod')
    if exist('cfg','var') & isfield(cfg,'MAIN') & isfield(cfg.MAIN,'auto') & cfg.MAIN.auto == 1
        settings.COV.covmethod = 'Folder';
    else
        disp ('   Ask for covariance-method')
        settings.COV = lab_get_cov([],header);
        if isempty(settings)
            return
        end
    end
end

if isfield(settings,'searchpath')
    searchpath = settings.searchpath;
    [~,~,searchformat,searchfileS] = lab_filename(settings.searchfile);
    settings = rmfield(settings,'searchpath');
    settings = rmfield(settings,'searchfile');
    clearvars tmp
elseif isfield(header,'EEG_filepath')
    searchpath = header.EEG_filepath;
    [~,~,searchformat,searchfileS] = lab_filename(header.EEG_file);
    clearvars tmp
else
    searchpath = [];
end

% Reduce to datachannels
if isfield(header,'numdatachannels')
    [data,header] = lab_reduce_channels(data,header,1:header.numdatachannels);
end

% Get Covariance from settings (COV-File)
if strcmp(settings.COV.covmethod,'COV-File') & isfield(settings.COV,'covfile') & ~isempty(settings.COV.covfile)
    if size(settings.COV.covfile,1) == size(data,1)
        disp ('   Use COV-File from settings')
        header.cov = settings.COV.covfile;
        return
    else
        disp ('   Warning: COV-File from settings not matching, use calculation from Input-File')
        settings.COV.covmethod = 'Input-File';
    end
end

% Collect data for covariance of all eeg/meg files in same folder
if ~isempty(searchpath) & strcmp(settings.COV.covmethod,'Folder')
    if exist(fullfile(searchpath,'covariance.cov'),'file')
        disp ('   Load covariance from previous run (covariance.cov)')
        header.cov = lab_read_cov(fullfile(searchpath,'covariance.cov'));
        return
    end
    files = dir(fullfile(searchpath,['*.' searchformat]));
    data2 = [];
    events.POS = [];
    events.DUR = [];
    events.TYP = [];
    for i = 1:size(files,1)
        headertmp = [];
        datatmp = [];
        try %#ok<TRYNC>
            [datatmp,headertmp,cfg] = lab_read_data(fullfile(searchpath,files(i,1).name),cfg);
        end
        if isfield(headertmp,'numdatachannels')
            disp('     do preprocessing for file in same folder')
            if isfield(cfg,'FILT') & ~isempty(cfg.FILT)
                cfgfilt = cfg.FILT;
                [datatmp,headertmp,cfg.FILT] = lab_filter(datatmp,headertmp,cfgfilt,'novrb');
            end
            if isfield(cfg,'RESAMPLE') & isfield(cfg.RESAMPLE,'new_samplingrate') & ~isempty(cfg.RESAMPLE.new_samplingrate)
                [datatmp,headertmp] = lab_resample_data(datatmp,headertmp,cfg);
            end
            if isfield(header,'includechans')
                [datatmp,headertmp] = lab_reduce_channels(datatmp,headertmp,header.includechans);
            end
            if isfield(header,'interpolated')
                if isfield(headertmp,'badchans')
                    headertmp.badchans = union(headertmp.badchans,header.interpolated);
                    header.badchans = header.badchans (:)';
                else
                    headertmp.badchans = header.interpolated;
                end
                headertmp.goodchans = setdiff(1:headertmp.numdatachannels,headertmp.badchans);
                headertmp.goodchans = headertmp.goodchans(:)';
                [datatmp,headertmp] = lab_interpolate_bad(datatmp,headertmp,'spherical');
            end
            datatmp = datatmp(1:headertmp.numdatachannels,:);
        end
        if size(data,1) == size(datatmp,1)
            if isfield(headertmp,'events') & ~isempty(header.events.POS)
                events.POS = [events.POS (headertmp.events.POS + size(data2,2))];
                events.DUR = [events.DUR headertmp.events.DUR];
                events.TYP = [events.TYP headertmp.events.TYP];
            end
            data2 = [data2 datatmp]; %#ok<AGROW>
        end
    end
    if ~isempty(data2)
        data = data2;
        header.events = events;
    end
    clearvars data2 events
end

% Limit data to predefined range
if isfield(settings,'range') & ~isempty(settings,'range')
    if settings.range(1) >= 1 & settings.range(1) < size(data,2)
        cfg.REDUCE.firstsample = settings.range(1);
    else
        cfg.REDUCE.firstsample = 1;
    end
    if length(settings.range) > 1 & settings.range(2) > cfg.REDUCE.firstsample & settings.range(2) <= size(data,2)
        cfg.REDUCE.lastsample = settings.range(2);
    else
        cfg.REDUCE.lastsample = size(data,2);
    end
else
    cfg.REDUCE.firstsample = 1;
    cfg.REDUCE.lastsample = size(data,2);
end
[data,header] = lab_reduce_datalength(data,header,cfg);

% Exclude markers
if isfield(header,'events') & ~isempty(header.events.POS) & ~isempty(settings.COV.markerexclude)
    events = header.events;
    for i = 1:size(settings.COV.markerexclude,1)
        markers = [];
        if strcmp(settings.COV.markerexclude{i,1},'all')
            tmp = 1:size(header.events.POS,2);
        else
            tmp = find(strcmp(events.TYP,settings.COV.markerexclude{i,1}));
        end
        markers = union(markers,tmp);
        clearvars tmp
    end
    for i = 1:length(markers)
        if events.POS(1,markers(i)) <= 1 & (events.POS(1,markers(i)) + events.DUR(1,markers(i))) > 1
            data = data(:,(events.POS(1,markers(i))+events.DUR(1,markers(i))):end);
        elseif (events.POS(1,markers(i))+events.DUR(1,markers(i))) >= size(data,2)
            data = data(:,1:(events.POS(1,markers(i))-1));
        elseif events.POS(1,markers(i)) > 1 & (events.POS(1,markers(i))+events.DUR(1,markers(i))) < size(data,2)
            data = [data(:,1:(events.POS(1,markers(i))-1)) data(:,(events.POS(1,markers(i))+events.DUR(1,markers(i))):end)];
        end
    end
end

% Include markers
if isfield(header,'events') & ~isempty(header.events.POS) & ...
        isfield(settings.COV,'markerinclude') & ~isempty(settings.COV.markerinclude)
    events = header.events;
    for i = 1:size(settings.COV.markerinclude,1)
        markers = [];
        if strcmp(settings.COV.markerinclude{i,1},'all')
            tmp = 1:size(header.events.POS,2);
        else
            tmp = find(strcmp(events.TYP,settings.COV.markerinclude{i,1}));
        end
        markers = union(markers,tmp);
        clearvars tmp
    end
    datatmp = [];
    for i = 1:length(markers)
        if events.POS(1,markers(i)) <= 1 & (events.POS(1,markers(i)) + events.DUR(1,markers(i))) > 1
            datatmp = [datatmp data(:,1:(events.POS(1,markers(i)) + events.DUR(1,markers(i))))]; %#ok<AGROW>
        elseif events.POS(1,markers(i)) < size(data,2) & (events.POS(1,markers(i)) + events.DUR(1,markers(i))) >= size(data,2)
            datatmp = [datatmp data(:,(events.POS(1,markers(i))):end)]; %#ok<AGROW>
        elseif events.POS(1,markers(i)) > 1 & (events.POS(1,markers(i))+events.DUR(1,markers(i))) < size(data,2)
            datatmp = [datatmp data(:,events.POS(1,markers(i)):(events.POS(1,markers(i)) + events.DUR(1,markers(i))))]; %#ok<AGROW>
        end
    end
    if ~isempty(datatmp)
        data = datatmp;
    else
        disp('     Warning: No markers to include, take whole input-file')
    end
end

% Calculate covariance matrix
if ~isempty(data)
    header.cov = cov(data(1:header.numdatachannels,:)');
end

% Write covariance
if ~isempty(searchpath)
    if strcmp(settings.COV.covmethod,'Folder')
        COV_file = fullfile(searchpath,'covariance.cov');
    else
        COV_file = fullfile(searchpath,[searchfileS '.cov']);
    end
    lab_write_cov(COV_file,header);
end

end