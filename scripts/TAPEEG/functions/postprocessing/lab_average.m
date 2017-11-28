% Average data on marker (eg for evoked potentials)
%
% [AVG,cfg] = lab_average(data,header,cfg)
%
% data                = matrix (channels x timeframes)
% header              = output of 'lab_read_data'
% cfg.Output_file     = filename of result
% cfg.Output_filepath = filepath of result
%
% Written by F. Hatz 2012 Neurology Basel

function [AVG,cfg] = lab_average(data,header,cfg)

AVG = [];
if ~exist('cfg','var') | ~isfield(cfg,'settings_path')
    cfg.settings_path = pwd;
end
if ~exist('header','var') | ~isfield(header,'samplingrate')
    header = lab_create_header(data);
end

if ~isfield(cfg,'Output_file')
    cfg.Output_file = header.EEG_file;
    cfg.Output_filepath = header.EEG_filepath;
    cfg.Output_fileS = lab_filename(cfg.Output_file);
end

if isfield(cfg,'SKIP') & isfield(cfg.SKIP,'AVG') & cfg.SKIP.AVG == true;
    return
end

disp('Averaging data')
if ~isfield(cfg,'AVG') | ~isfield(cfg.AVG,'marker')
    [cfg,skipprocessing] = lab_set_average(cfg,header,data);
    if skipprocessing == 1
        return
    else
        pause(0.1);
    end
end

% look for markers
if strcmp(cfg.AVG.marker{1,1},'*') | isempty(cfg.AVG.marker)
    if isfield(cfg,'MAIN') & isfield(cfg.MAIN,'auto') & cfg.MAIN.auto == 1
        if ~exist(fullfile(cfg.settings_path,'marker.xls'),'file')
            markerselection = unique(header.events.TYP);
            lab_write_xls(fullfile(cfg.settings_path,'marker~.xls'),markerselection');
            clearvars markerselection
        end
        while ~exist(fullfile(cfg.settings_path,'marker.xls'),'file')
            disp('  Pease edit -marker~.xls- and rename to -marker.xls-')
            pause(30);
        end
        if ispc
            [~,cfg.AVG.marker] = xlsread(fullfile(cfg.settings_path,'marker.xls'));
        else
            [~,cfg.AVG.marker] = xlsread(fullfile(cfg.settings_path,'marker.xls'),1,'','basic');
        end
        clearvars tmp
    else
        disp ('Ask for Markers')
        strlist = unique(header.events.TYP);
        defaultanswer = [];
        for i = 1:size(cfg.AVG.marker,1)
            defaultanswer = [defaultanswer find(strcmp(strlist,cfg.AVG.marker(i,1)))]; %#ok<AGROW>
        end
        [selection] = listdlg('PromptString','Select Markers:','SelectionMode','multiple', ...
            'ListString',strlist,'InitialValue',defaultanswer);
        cfg.AVG.marker = strlist(selection);
        clearvars selection strlist
        pause(0.2);
    end
end

AVGmarker = cfg.AVG.marker;
for i = 1:length(AVGmarker)
    tmp = header.events.POS(1,ismember(header.events.TYP,AVGmarker{i})==1);
    tmp = tmp - int64(cfg.AVG.markerOffset);
    tmp = tmp(1,tmp > 0);
    markerevents{i} = tmp(1,tmp < (size(data,2) - cfg.AVG.markerlength)); %#ok<AGROW>
end
if isfield(cfg.AVG,'combinemarker') & cfg.AVG.combinemarker == true & length(AVGmarker) > 1
    tmp = [];
    for i = 1:length(markerevents)
        tmp = union(tmp,markerevents{i});
    end
    if size(tmp,1) > 1
        tmp = tmp';
    end
    clearvars markerevents
    markerevents{1} = tmp;
    AVGmarker = cellstr(sprintf('%s_',cfg.AVG.marker{:}));
end

% Reduce data to 'cfg.AVG.includechans'
if isfield(cfg.AVG,'includechans') & ~isempty(cfg.AVG.includechans)
    [data,header] = lab_reduce_channels(data,header,cfg.AVG.includechans);
end

% Filter data
if isfield(cfg.AVG,'FILT') & isfield(cfg.AVG.FILT,'filtermode')
    [data,header,cfg.AVG.FILT] = lab_filter(data,header,cfg.AVG.FILT,'novrb');
end

% Correct montage
if strcmp(cfg.AVG.eegsource,'montage')
    if ~isfield(cfg.AVG,'montage') | isempty(cfg.AVG.montage)
        disp('   no montage loaded, use input structure')
        cfg.AVG.montage = lab_create_montage(size(data,1),header);
    end
    if cfg.AVG.montage(1,1).numchans ~= header.numchannels
        cfg.AVG.montage = lab_reduce_montage(cfg.AVG.montage,cfg,header,true);
    end
    if cfg.AVG.montage(1,1).numchans > header.numchannels
        disp('   invalid montage loaded, use input structure')
        montage = lab_create_montage(size(data,1),header);
    else
        montage = cfg.AVG.montage;
    end
end

% Average data
if size(markerevents{1},2) == 0
    disp('   Average not possible')
    return
end

% Create Output folder
if ~isfield(cfg.AVG,'folder')
    cfg.AVG.folder = 'AVG';
end
warning off %#ok<WNOFF>
mkdir(fullfile(cfg.Output_filepath,cfg.AVG.folder));
warning on %#ok<WNON>
Average_filepath = fullfile(cfg.Output_filepath,cfg.AVG.folder);
Output_fileS = cfg.Output_fileS;

cfg.SKIP.AVG = true;
SKIP = cfg.SKIP;
for j = 1:length(markerevents)
    if isnumeric(cfg.AVG.eegsource) | ~isnan(str2double(cfg.AVG.eegsource))
        AVGsource = 'chanref';
    else
        AVGsource = cfg.AVG.eegsource;
    end
    if strcmp(AVGsource,'montage')
        if cfg.AVG.interpolate == 1
            disp('   Interpolate bad channels')
            [datacalc,headercalc] = lab_interpolate_bad(data,header);
        end
        [datacalc,headercalc] = lab_references(datacalc,headercalc,montage,cfg.AVG);
    elseif strcmp(AVGsource,'mean') | strcmp(AVGsource,'median') | strcmp(AVGsource,'laplacian')
        [datacalc,headercalc] = lab_references(data,header,AVGsource,cfg.AVG);
        if cfg.AVG.interpolate == 1
            disp('   Interpolate bad channels')
            [datacalc,headercalc] = lab_interpolate_bad(datacalc,headercalc);
        end
    elseif strcmp(AVGsource,'chanref')
        [datacalc,headercalc] = lab_references(data,header,cfg.AVG.eegsource,cfg.AVG);
        if cfg.AVG.interpolate == 1
            disp('   Interpolate bad channels')
            [datacalc,headercalc] = lab_interpolate_bad(datacalc,headercalc);
        end
    else
        datacalc = data;
        headercalc = header;
        if cfg.AVG.interpolate == 1
            disp('   Interpolate bad channels')
            [datacalc,headercalc] = lab_interpolate_bad(datacalc,headercalc);
        end
    end
    if isempty(datacalc)
        return
    end
    
    if length(markerevents) > 1
        Average_file = [Output_fileS '_T' num2str(j)];
    else
        Average_file = Output_fileS;
    end
    
    [average,epochs,epochsValid,epochsBad] = lab_average_epochs(datacalc,headercalc,markerevents{j},cfg.AVG);
    if max(average(:)) == 0
        disp('     Abort avergaing, no valid sweeps found')
        return
    end
    
    if isfield(cfg.AVG,'JITTER') & ~isempty(cfg.AVG.JITTER)
        Jitter = lab_calculate_jitter(average,epochs,epochsValid,cfg.AVG.JITTER);
    end
    if isfield(cfg.AVG,'Correct') & isfield(cfg.AVG.Correct,'BAD')
        disp('   Detect bad channels (average)')
        headercalc.badchans = lab_detect_bad(average,headercalc,cfg.AVG.Correct.BAD,'noverbose');
    end
    if isfield(cfg.AVG,'Correct') & isfield(cfg.AVG.Correct,'method') & ~strcmp(AVGsource,'montage')
        disp('   Interpolation of bad channels (average)')
        [average,headercalc] = lab_interpolate_bad(average,headercalc,cfg.AVG.Correct.method);
    end
    electrodesvalid = sum(epochsValid,2) ./ size(epochsValid,2);
    channelstmp = cellstr(headercalc.channels);
    electrodesvalid = [channelstmp num2cell(electrodesvalid)];
    electrodesvalid = cat(1,[{'channel'} {'% good epochs'}],electrodesvalid);
    clearvars channelstmp
    
    % write result
    disp (['   Write result ' Average_file])
    Average_file = [Average_file '_Avg']; %#ok<AGROW>
    headercalc.numtimeframes = size(average,2);
    headercalc.events.POS = int64(1 + cfg.AVG.markerOffset);
    headercalc.events.DUR = int64(0);
    headercalc.events.OFF = int64(0);
    headercalc.events.TYP = {'Trigger'};
    AVG = cfg.AVG;
    save(fullfile(Average_filepath,[Average_file '.mat']),'epochsValid','AVG','electrodesvalid','epochsBad');
    clearvars AVG
    Filename = fullfile(Average_filepath,[Average_file '.sef']);
    for nfor = 1:length(cfg.AVG.format)
        [~,cfg.AVG] = lab_save_data(average,headercalc,cfg.AVG.format{nfor},Filename,cfg.AVG,cfg);
    end
    lab_write_xls([Filename(1:end-4) '.xls'],electrodesvalid);
    if exist('Jitter','var') & ~isempty(Jitter)
        xlsout = [cellstr(header.channels(1:size(average,1),:)) num2cell(Jitter.ChanValid)' num2cell(Jitter.Mean) num2cell(Jitter.Std) num2cell(Jitter.MeanError)];
        xlsout = cat(1,{'Channel','Valid','Lag_Mean','Lag_Std','Error_Mean'},xlsout);
        lab_write_xls([Filename(1:end-4) '_Jitter.xlsx'],xlsout);
        clearvars xlsout
        xlsout = cellstr('Channel');
        for xj = 1:size(Jitter.Lag,2)
            xlsout{1,xj+1} = ['Sweep' num2str(xj)];
        end
        xlsout = cat(1,xlsout,[cellstr(header.channels(1:size(average,1),:)) num2cell(Jitter.Lag)]);
        lab_write_xls([Filename(1:end-4) '_JitterLags.xlsx'],xlsout);
        clearvars xlsout xj
    end
    epochsValid = mean(epochsValid,1);
    if length(unique(epochsValid)) == 2 & min(epochsValid) == 0 & max(epochsValid) == 1
        events.POS = markerevents{j}(epochsValid == 1);
        events.TYP = repmat({'Sweep'},1,length(events.POS));
        events.DUR = repmat(int64(cfg.AVG.markerlength),1,length(events.POS));
        events.OFF = repmat(int64(0),1,length(events.POS));
        lab_write_mrk([Filename(1:end-4) 'Select.mrk'],events);
        clearvars events
    end
    clearvars epochsValid headerAvg electrodesvalid Filename
    
    % Store result in Output variable
    AVG.(AVGmarker{j}) = average;
    if exist('Jitter','var') & ~isempty(Jitter)
        AVG.([AVGmarker{j} '_Jitter']) = Jitter;
    end
    
    % Dipol Fit
    if isfield(cfg.AVG,'IS') & ~isempty(cfg.AVG.IS) & ...
            (strcmp(cfg.AVG.IS.eegsource,AVGsource) | strcmp(cfg.AVG.IS.eegsource,'input'))
        warning off %#ok<WNOFF>
        mkdir(fullfile(Average_filepath,'DipolFit'));
        cfg.AVG.Output_filepath = fullfile(Average_filepath,'DipolFit');
        cfg.AVG.Output_file = [Average_file(1:end-4) '.sef'];
        warning on %#ok<WNON>
        if isfield(cfg.AVG.IS,'averageepochs') & cfg.AVG.IS.averageepochs == 1
            [~,cfg.AVG] = lab_inversesolution(average,headercalc,cfg.AVG);
        else
            [~,cfg.AVG] = lab_inversesolution(epochs,headercalc,cfg.AVG);
        end
    end
    
    % -- Do calculations on average --
    Output_file = cfg.Output_file;
    Output_fileS = cfg.Output_fileS;
    cfg.Output_file = [Average_file '.sef'];
    cfg.SKIP = SKIP;
    cfg = lab_postprocessing(average,headercalc,cfg,true);
    cfg.Output_file = Output_file;
    cfg.Output_fileS = Output_fileS;
    
    % write verbose
    Verbose_file=[Average_file '.vrb'];
    fid=fopen(fullfile(Average_filepath,Verbose_file),'w');
    fprintf(fid,'Average\n');
    fprintf(fid,datestr(now,0));
    fprintf(fid,'\n');
    fprintf(fid,'Files:\n');
    fprintf(fid,['EEG input file: ' cfg.Output_file]);
    fprintf(fid,'\n');
    if strcmp(AVGsource,'laplacian')
        fprintf(fid,'Laplacian Maxdistance: ');
        fprintf(fid,num2str(cfg.REF.lap_maxdistance));
        fprintf(fid,'\n');
        fprintf(fid,'Laplacian Weight Maxdistance: ');
        fprintf(fid,num2str(cfg.REF.lap_weightmaxdistance));
        fprintf(fid,'\n');
    end
    if isfield(cfg.AVG,'FILT') & ~isempty(cfg.AVG.FILT)
        fprintf(fid,'Filter settings: ');
        fprintf(fid,cfg.AVG.FILT.filtermode);
        fprintf(fid,'\n');
        fprintf(fid,'  Notch Filter: ');
        for f = 1:size(cfg.AVG.FILT.notch,2)
            fprintf(fid,[num2str(cfg.AVG.FILT.notch(1,f)) 'Hz ']);
        end
        fprintf(fid,'\n');
        fprintf(fid,['  Hghpass Filter: ' num2str(cfg.AVG.FILT.highpass) 'Hz']);
        fprintf(fid,'\n');
        fprintf(fid,['  Lowpass Filter: ' num2str(cfg.AVG.FILT.lowpass) 'Hz']);
        fprintf(fid,'\n\n');
    end
    if cfg.AVG.interpolate == 1
        fprintf(fid,'Interpolate bad channels: enabled');
    else
        fprintf(fid,'Interpolate bad channels: disabled');
    end
    fprintf(fid,'\n');
    fprintf(fid,['Marker: ' AVGmarker{j}]);
    fprintf(fid,'\n');
    fprintf(fid,'Length of epochs: ');
    fprintf(fid,num2str(cfg.AVG.markerlength));
    fprintf(fid,'\n');
    fprintf(fid,'Length of pre-run: ');
    fprintf(fid,num2str(cfg.AVG.markerOffset));
    fprintf(fid,'\n');
    fprintf(fid,['Method for average: ' cfg.AVG.AVGmethod]);
    fprintf(fid,'\n');
    if isfield(cfg.AVG.Reject,'BAD')
        if isnumeric(cfg.AVG.Reject.BAD)
            fprintf(fid,['Standard deviation for reject: ' num2str(cfg.AVG.Reject.BAD)]);
        elseif isstruct(cfg.AVG.Reject.BAD)
            fprintf(fid,'Settings for detecting bad epochs:\n');
            lab_write_bad_vrb(fid,cfg.AVG.Reject.BAD,epochsBad);
            fprintf(fid,'\n\n');
        end
        if strcmp(cfg.AVG.Reject.method,'single')
            fprintf(fid,' (single channels)');
        elseif strcmp(cfg.AVG.Reject.method,'percent')
            fprintf(fid,[' (' num2str(cfg.AVG.Reject.percent) '% of channels)']);
        end
    else
        fprintf(fid,'Standard deviation for reject: disabled');
    end
    fprintf(fid,'\n');
    if isfield(cfg.AVG,'Correct') & isfield(cfg.AVG.Correct,'BAD') & ~isempty(cfg.AVG.Correct.BAD)
        fprintf(fid,'Bad channels detection after averaging: enabled');
    else
        fprintf(fid,'Bad channels detection after averaging: disabled');
    end
    fprintf(fid,'\n');
    if isfield(cfg.AVG,'Correct') & isfield(cfg.AVG.Correct,'method') & ...
            ~isempty(cfg.AVG.Correct.method) & ~strcmp(cfg.AVG.Correct.method,'disabled')
        fprintf(fid,'Bad channels interpolation after averaging: enabled');
    else
        fprintf(fid,'Bad channels interpolation after averaging: disabled');
    end
    fprintf(fid,'\n');
    fclose(fid);
end

return