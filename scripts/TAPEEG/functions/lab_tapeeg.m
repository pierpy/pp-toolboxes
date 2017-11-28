% TAPEEG: Main_App
%
% Written by F. Hatz, Neurology University Basel 2012
%
% cfg.SEARCH.searchfolder   =   folder to process
% cfg.SEARCH.searchstring   =   searchstring (eg '.edf')
% cfg.MAIN.log_path         =   path to log-file
%
% to start ICA-backtransform select *_ICA.mat for file or search string

function cfg = lab_tapeeg(cfg)

global VERSION Settings_Path Segments
if isempty(Settings_Path)
    Settings_Path = pwd;
end

calc.Filelist = [];

disp('--------------')
disp(['TAPEEG ' VERSION])
disp('--------------')

if exist('cfg','var') & isfield(cfg,'MAIN')
    cfg.MAIN.auto = 1;
else
    cfg.MAIN.auto = 0;
end

% Empty all global variables 
lab_empty_global_vars;

% select single file or folder
if cfg.MAIN.auto == 0 | ~isfield(cfg.SEARCH,'searchfolder') | isempty(cfg.SEARCH.searchfolder) | ...
        ~isfield(cfg.SEARCH,'searchstring') | isempty(cfg.SEARCH.searchstring)
    disp('Select mode')
    button = questdlg('Single File or Folder','Select mode','Cancel','Folder','File','File');
    if strcmp(button,'Folder')
        disp ('Ask for folder')
        cfg.SEARCH.searchfolder = uigetdir(pwd,'Select folder');
        if isempty(cfg.SEARCH.searchfolder) | cfg.SEARCH.searchfolder == 0
            return
        else
            cd(cfg.SEARCH.searchfolder);
        end
    elseif strcmp(button,'File')
        disp ('Ask for EEG-file')
        [EEG_file,EEG_filepath]=uigetfile('*.*','Select EEG file');
        if ~isempty(EEG_file) & EEG_file ~= 0
            cd(EEG_filepath);
            if exist(fullfile(EEG_filepath,'signal1.bin'),'file')
                tmp = strfind(EEG_filepath,filesep);
                EEG_file = [EEG_filepath(tmp(end-1)+1:tmp(end)-1) '.mff'];
                EEG_filepath = EEG_filepath(1:tmp(end-1));
                clearvars tmp
            elseif exist(fullfile(fullfile(EEG_filepath,EEG_file),'signal1.bin'),'file')
                EEG_file = [EEG_file '.mff'];
            end
            calc.Filelist{1,1} = fullfile(EEG_filepath,EEG_file);
            calc.Filelist_done = [];
            calc.Filelist_doneall = [];
            cfg.SEARCH.searchfolder = 0;
            clearvars EEG_file EEG_filepath
        else
            return
        end
    else
        return
    end
end

% set settings path
if isfield(cfg.SEARCH,'searchfolder') & ischar(cfg.SEARCH.searchfolder) & exist(cfg.SEARCH.searchfolder,'dir')
    if strcmp(cfg.SEARCH.searchfolder(end),filesep)
        cfg.SEARCH.searchfolder = cfg.SEARCH.searchfolder(1:end-1);
    end
    if ~isfield(cfg,'SettingsFileAuto')
        cfg.settings_path = cfg.SEARCH.searchfolder;
    end
elseif ~isempty(calc.Filelist)
    pathtmp = calc.Filelist{1,1};
    tmp=strfind(pathtmp,filesep);
    if size(tmp,2) > 0 & exist(fullfile(pathtmp(1:tmp(end)),'settings.mat'),'file')
        cfg.settings_path = pathtmp(1:tmp(end));
    elseif size(tmp,2) > 1 & exist(fullfile(pathtmp(1:tmp(end-1)),'settings.mat'),'file')
        cfg.settings_path = pathtmp(1:tmp(end-1));
    elseif size(tmp,2) > 2 & exist(fullfile(pathtmp(1:tmp(end-2)),'settings.mat'),'file')
        cfg.settings_path = pathtmp(1:tmp(end-2));
    else
        cfg.settings_path = pathtmp(1:tmp(end));
        if exist(fullfile(cfg.settings_path,'signal1.bin'),'file')
            cfg.settings_path = pathtmp(1:tmp(end-1));
        end
    end
    clearvars pathtmp tmp
else
    cfg.settings_path = pwd;
end

% set name of settings_file
if ~isfield(cfg,'settings_file')
    cfg.settings_file = 'settings.mat';
end

% set log-file
if ~isfield(cfg.MAIN,'log_path')
    cfg.MAIN.log_path = fullfile(cfg.settings_path,'TAPEEG.log');
end

if cfg.MAIN.auto == 0
    if exist(cfg.MAIN.log_path,'file')
        [Log_file,Log_filepath] = lab_filename(cfg.MAIN.log_path);
        Settings_path = fullfile(Log_filepath,'settings');
        if ~exist(Settings_path,'dir')
            mkdir(Settings_path);
        end
        if exist(fullfile(Log_filepath,'settings.mat'),'file')
            tmp = dir(fullfile(Log_filepath,'settings.mat'));
        else
            tmp = dir(cfg.MAIN.log_path);
        end
        tmp = textscan(tmp.date,'%s ');
        Log_date = tmp{1,1}{1,1};
        clearvars tmp
        Settings2_path = fullfile(Settings_path,Log_date);
        j = 1;
        while exist(Settings2_path,'dir')
            Settings2_path = fullfile(Settings_path,[Log_date '_' num2str(j)]);
            j = j + 1;
        end
        clearvars j
        
        warning off %#ok<WNOFF>
        mkdir(Settings2_path);
        warning on %#ok<WNON>
        try %#ok<TRYNC>
            movefile(cfg.MAIN.log_path,fullfile(Settings2_path,Log_file));
        end
        try %#ok<TRYNC>
            if exist(fullfile(cfg.settings_path,'settings.mat'),'file')
                movefile(fullfile(cfg.settings_path,'settings.mat'),fullfile(Settings2_path,'settings.mat'));
                cfg.settings_fileM = fullfile(Settings2_path,'settings.mat');
            elseif exist(fullfile(Log_filepath,'settings.mat'),'file')
                movefile(fullfile(Log_filepath,'settings.mat'),fullfile(Settings2_path,'settings.mat'));
                cfg.settings_fileM = fullfile(Settings2_path,'settings.mat');
            end
        end
        pause(0.5);
    end
    clearvars Log_file Log_filepath Log_date Settings_path Settings2_path
end

% start Log-file
diary(cfg.MAIN.log_path);

% open settings-file
if ~isfield(cfg,'SettingsFileAuto')
    try
        warning off %#ok<WNOFF>
        [cfg,skipprocessing] = lab_set_settings(cfg);
        if skipprocessing == 1;
            return
        end
        warning on %#ok<WNON>
    catch %#ok<CTCH>
        disp ('Loading settings-file failed')
    end
    cfg.SettingsFileAuto = false;
    if ~isfield(cfg,'settings_file') | isempty(cfg.settings_file)
        cfg.settings_file = 'settings.mat';
    end
end
Settings_Path = cfg.settings_path;

% Search files if needed
if isfield(cfg.SEARCH,'searchfolder') & ~cfg.SEARCH.searchfolder == 0
    if cfg.MAIN.auto == 1
        if ~isfield(cfg.SEARCH,'searchstring') | isempty(cfg.SEARCH.searchstring)
            disp('Abort, no valid search-string set')
            return
        end
    else
        [cfg,skipprocessing] = lab_set_searchstrings(cfg,0,1);
        if skipprocessing == 1
            diary off
            return
        end
    end
    if ~isempty(cfg.SEARCH.searchstring)
        [calc,cfg] = lab_search_files(cfg);
    else
        calc.Filelist = [];
    end
end
if size(calc.Filelist,2) == 0
    disp('No files to process')
    diary off
    return
end

if cfg.MAIN.auto == 1 | (isfield(cfg,'editsettings') & cfg.editsettings == false)
    skipsettings = 1;
else
    skipsettings = 0;
end

% set flag for last files in folders
[calc.Filefirst,tmp] = lab_find_lastperfolder(calc.Filelist);
calc.Filelast = tmp;
clearvars tmp;

for filenr = 1:size(calc.Filelist,2)
    % Process files
    if isfield(cfg,'patient')
        cfg = rmfield(cfg,'patient');
    end
    [cfg.EEG_file,cfg.EEG_filepath,cfg.fileformat,cfg.EEG_fileS] = lab_filename(calc.Filelist{1,filenr});
    if length(cfg.EEG_file) > 8 & strcmp(cfg.EEG_file(end-7:end),'_ICA.mat')
        doICAback = true;
    else
        doICAback = false;
    end
    if strcmp(cfg.EEG_fileS(end-3:end),'filt')
        skipfilt = 1;
    else
        skipfilt = 0;
    end
    if filenr == size(calc.Filelist,2)
        cfg.lastfile = true;
    else
        cfg.lastfile = false;
    end
    if filenr == 1
        cfg.firstfile = true;
    else
        cfg.firstfile = false;
    end
    cfg.lastfilefolder = calc.Filelast(filenr);
    cfg.firstfilefolder = calc.Filefirst(filenr);
    
    skipprocessing = 0;
    skiperror = 0; %#ok<NASGU>
    skippreprocessing = 0;
    skippostprocessing = 0;
    cd(cfg.EEG_filepath);
    
    if doICAback == true
        % Start ICA backtransform
        if ~isfield(cfg,'ICABACK') | isempty(cfg.ICABACK) | skipsettings == 0
            header = lab_load_ICAheader(fullfile(cfg.EEG_filepath,cfg.EEG_file));
            [cfg,~,skipprocessing] = lab_edit_cfg(cfg,header,1,1);
            if skipprocessing == 1
                diary off
                return
            else
                skipsettings = 1;
            end
        end
        if isfield(cfg,'ICABACK') & ~isempty(cfg.ICABACK)
            [data,header,cfg,skipprocessing] = lab_ICAback(fullfile(cfg.EEG_filepath,cfg.EEG_file),cfg);
            if isfield(header,'patient') & cfg.MAIN.auto == 1
                cfg.patient = header.patient;
            else
                [~,cfg,skipprocessing] = lab_subjectname(fullfile(cfg.EEG_filepath,cfg.EEG_file),cfg);
                if skipprocessing == 1
                    diary off
                    return
                end
                header.patient = cfg.patient;
            end
        else
            data = [];
            header = [];
        end
        if ~isempty(data)
            [~,~,cfg.fileformat,cfg.EEG_fileS] = lab_filename(cfg.EEG_file);
            cfg.SEG.select = 'Select complete file';
            skippreprocessing = 1;
        else
            skipprocessing = 1;
            skiperror = 1; %#ok<NASGU>
        end
    else
        % read eeg file
        disp(['Read file (' cfg.EEG_file ')'])
        if cfg.MAIN.auto == 1 & cfg.SettingsFileAuto == true
            disp('    Control if file is still created')
            controlwriting = 1;
            bytes = [];
            while controlwriting == 1
                if strcmp(cfg.fileformat,'mff')
                    tmp = dir(fullfile(fullfile(cfg.EEG_filepath,cfg.EEG_fileS),'signal1.bin'));
                else
                    tmp = dir(fullfile(cfg.EEG_filepath,cfg.EEG_file));
                end
                if tmp.bytes == bytes
                    controlwriting = 0;
                end
                bytes = tmp.bytes;
                pause(4);
            end
            clearvars controlwriting bytes
        end
        
        [data,header,cfg] = lab_read_data(fullfile(cfg.EEG_filepath,cfg.EEG_file),cfg,true);
        if strcmp(header.datatype,'matrix')
            if isempty(data)
                [data,header,cfg] = lab_read_data(fullfile(cfg.EEG_filepath,cfg.EEG_file),cfg,false,true);
            end
            if filenr == size(calc.Filelist,2)
                cfg.lastfile = true;
            else
                cfg.lastfile = false;
            end
            if skipsettings == 0
                cfg.Output_filepath = cfg.EEG_filepath;
                cfg.Output_file = cfg.EEG_file;
                if size(data,3) > 1
                    if isfield(header,'subjects') & ~isempty(header.subjects) & ...
                            iscell(header.subjects) & ~isempty(header.subjects{1});
                        cfg.Output_file = header.subjects{1};
                    end
                    [cfg,skipprocessing] = lab_set_process_matrix(cfg,false,false,false,true,header);
                else
                    [cfg,skipprocessing] = lab_set_process_matrix(cfg,false,false,false,false,header);
                end
                if skipprocessing == 1
                    return
                else
                    skipsettings = 1;
                end
            end
            if isfield(cfg,'MATRIX') & ~isempty(cfg.MATRIX)
                cfg.Output_filepath = cfg.EEG_filepath;
                cfg.Output_file = cfg.EEG_file;
                [~,~,cfg] = lab_process_matrix(data,header,cfg);
            else
                disp('  skip processing matrix, no settings')
            end
            skipprocessing = 1;
            if ~exist('flagsettings','var') & cfg.SettingsFileAuto == false
                lab_save_settings(cfg,calc);
                flagsettings = 1; %#ok<NASGU>
            end
        elseif strcmp(header.datatype,'mri')
            if skipsettings == 0
                [cfg,skipprocessing] = lab_set_process_mri(cfg);
                if skipprocessing == 1
                    return
                else
                    skipsettings = 1;
                end
            end
            if isfield(cfg,'MRI') & ~isempty(cfg.MRI)
                cfg = lab_process_mri(cfg);
            else
                disp('  skip processing mri, no settings')
            end
            skipprocessing = 1;
            if ~exist('flagsettings','var') & cfg.SettingsFileAuto == false
                lab_save_settings(cfg,calc);
                flagsettings = 1; %#ok<NASGU>
            end
        elseif strcmp(header.datatype,'eeg')
            if ~isfield(header,'numdatachannels')
                header.numdatachannels = header.numchannels;
            end
            
            % check for constant number of channels
            if exist('settingschannels','var') & cfg.MAIN.auto == 0
                if settingschannels ~= header.numdatachannels
                    if ~isfield(cfg,'LOCS') | ~isfield(cfg.LOCS,'matchlocs') | cfg.LOCS.matchlocs == false
                        skipsettings = 0;
                    end
                end
            end
            
             % set settings (only first file)
            if ~isfield(cfg,'EXTRA') | skipsettings == 0
                if skipfilt == 1
                    [cfg,header,skipprocessing] = lab_edit_cfg(cfg,fullfile(cfg.EEG_filepath,cfg.EEG_file),1,0,1);
                else
                    [cfg,header,skipprocessing] = lab_edit_cfg(cfg,fullfile(cfg.EEG_filepath,cfg.EEG_file),1,0,0);
                end
                pause(0.1);
                if skipprocessing == 1
                    diary off
                    return
                elseif isfield(cfg,'ICA') & isfield(cfg.ICA,'automated') & ...
                        cfg.ICA.automated == false
                    skippostprocessing = 1;
                end
                skipsettings = 1;
                settingschannels = header.numdatachannels;
            end
            
            if isfield(cfg,'LOCS') & isfield(cfg.LOCS,'matchlocs') & ...
                    cfg.LOCS.matchlocs == true & ~exist('FlagMatchLocs','var')
                cfg = lab_find_matchchannels(cfg,calc);
                FlagMatchLocs = true; %#ok<NASGU>
                [~,header,cfg] = lab_read_data(fullfile(cfg.EEG_filepath,cfg.EEG_file),cfg,true);
            end
                        
            if ~isfield(cfg,'EXTRA')
                disp('Abort: wrong settings selected (extra channels info missing)')
                diary off
                return
            end
            
            [~,cfg,skipprocessing] = lab_subjectname(fullfile(cfg.EEG_filepath,cfg.EEG_file),cfg);
            if skipprocessing == 1
                diary off
                return
            end
            [~,~,~,cfg.EEG_fileS] = lab_filename(cfg.EEG_file);
            cfg.Output_file = [cfg.EEG_fileS '.sef'];
            if strcmp(cfg.fileformat,'mff') & length(cfg.EEG_filepath) < length(cfg.settings_path)
                if cfg.subjectname == -1
                    cfg.patient = cfg.EEG_fileS;
                end
                cfg.EEG_filepath = [fullfile(cfg.EEG_filepath,cfg.EEG_fileS) filesep];
            end
            cfg.Output_filepath = cfg.EEG_filepath;
        elseif strcmp(header.datatype,'error')
            skipprocessing = 1;
            skiperror = 1; %#ok<NASGU>
        end
    end
    
    if skipprocessing == 0
        % Define file reading
        if isfield(cfg,'FILEREADING') & ~isempty(cfg.FILEREADING) & doICAback == false
            if isempty(data)
                [Parts,cfg] = lab_define_filereading(header,cfg);
                if isempty(Parts)
                    skipprocessing = 1;
                end
            else
                disp('Warning: Segment-Reading not possible for this filetype')
            end
        else
            Parts = [1 header.numtimeframes];
        end
        
        EEG_file = cfg.EEG_file;
        EEG_filepath = cfg.EEG_filepath;
        Seg_Number = 1;
        for S = 1:size(Parts,1)
            cfg.EEG_file = EEG_file;
            cfg.EEG_filepath = EEG_filepath;
            if skipprocessing == 0
                % Define segments in eeg to process
                if doICAback == true
                    Segments{1,1} = data;
                    Segments{2,1} = header;
                elseif strcmp(cfg.EEG_fileS(end-3:end),'filt')
                    if isempty(data) | size(data,2) < 2
                        [data,header,cfg] = lab_read_data(fullfile(cfg.EEG_filepath,cfg.EEG_file),cfg,false,true);
                    end
                    skipfilt = 1;
                    Segments{1,1} = data;
                    Segments{2,1} = header;
                else
                    % Read EEG/MEG
                    if isempty(data) | size(data,2) < 2
                        [data,header,cfg] = lab_read_data(fullfile(cfg.EEG_filepath,cfg.EEG_file),cfg,false,Parts(S,:),true);
                    end
                    
                    % Edit Markers
                    if isfield(cfg,'MARK') & ~isempty(cfg.MARK)
                        [data,header,cfg,skipprocessing] = lab_edit_markers(data,header,cfg);
                        if skipprocessing == 1
                            diary off
                            return
                        end
                    end
                    
                    % Create Segments
                    if skipprocessing == 0
                        if ~exist('cfg','var') | ~isfield(cfg,'SEG') | isempty(cfg.SEG)
                            cfg.SEG.select = {'Select complete file'};
                        end
                        [Segments,cfg,skipprocessing] = lab_define_segments(data,header,cfg);
                        if skipprocessing == 1
                            diary off
                            return
                        end
                        if isempty(Segments)
                            lab_save_settings(cfg,calc);
                            skipprocessing = 1;
                            skiperror = 1; %#ok<NASGU>
                        end
                    end
                end
                header.patient = cfg.patient;
            end
            
            if skipprocessing == 0
                EEG_fileS = cfg.EEG_fileS;
                numsegments = size(Segments,2);
                for k = 1:numsegments
                    skipprocessing = 0;
                    skiperror = 0;
                    
                    % set flag for last segments
                    if filenr == size(calc.Filelist,2)
                        if k == numsegments
                            cfg.lastsegment = true;
                        else
                            cfg.lastsegment = false;
                        end
                    else
                        [~,TMP_filepath] = lab_filename(calc.Filelist{1,filenr});
                        [~,TMP_filepath2] = lab_filename(calc.Filelist{1,filenr+1});
                        if k == numsegments & ~strcmp(TMP_filepath,TMP_filepath2)
                            cfg.lastsegment = true;
                        else
                            cfg.lastsegment = false;
                        end
                    end
                    
                    % Process segments of eeg
                    data = Segments{1,k};
                    header = Segments{2,k};
                    header.numsegment = k;
                    if k == size(Segments,2)
                        StoreAll = 1;
                    else
                        StoreAll = 0;
                    end
                    
                    if size(Parts,1) > 1 | size(Segments,2) > 1
                        cfg.EEG_file = [EEG_fileS '_S' num2str(Seg_Number) '.' cfg.fileformat];
                        Seg_Number = Seg_Number + 1;
                        [~,~,~,cfg.EEG_fileS] = lab_filename(cfg.EEG_file);
                    else
                        clearvars Segments
                    end
                    cfg.Output_file = cfg.EEG_file;
                    cfg.Output_fileS = cfg.EEG_fileS;
                    
                    if ~isempty(data)
                        % Define name of processed eeg file
                        WorkingEEG_file = fullfile(cfg.EEG_filepath,cfg.EEG_file);
                        if skipfilt == 1
                            cfg.EEG_fileS = cfg.EEG_fileS(1:end-4);
                            [~,~,format] = lab_filename(cfg.EEG_file);
                            cfg.EEG_file = [cfg.EEG_fileS '.' format];
                            clearvars format
                            cfg.Output_file = [cfg.EEG_fileS '.sef'];
                        end
                        if max(strcmp(calc.Filelist_done,WorkingEEG_file)) & cfg.SEARCH.excludeprocessed == true
                            disp(['Skip ' cfg.EEG_file '(already processed)'])
                            skipprocessing = 1;
                        else
                            if numsegments > 1
                                disp(['Processing ' cfg.EEG_file '(' num2str(k) ' of ' num2str(numsegments) ')'])
                            else
                                disp(['Processing ' cfg.EEG_file])
                            end
                        end
                    else
                        disp('Abort: no data')
                        skipprocessing = 1;
                        skiperror = 1;
                    end
                    
                    if skipprocessing == 0
                        % remove locs if not matching
                        if isfield(header,'locs') & isfield(header.locs,'x') & length(header.locs.x) ~= header.numdatachannels
                            header = rmfield(header,'locs');
                            disp('   Delete locs, number of channels not matching')
                        end
                        
                        % control channels to exclude in results
                        if isfield(cfg,'numdatachans') & ~isempty(cfg.numdatachans) & cfg.numdatachans ~= header.numdatachannels
                            cfg = rmfield(cfg,'exclude');
                        end
                        if isfield(cfg,'exclude') & max(cfg.exclude) > header.numchannels
                            cfg = rmfield(cfg,'exclude');
                        end
                        if ~isfield(cfg,'exclude')
                            if cfg.MAIN.auto == 0
                                cfg = lab_set_exclude(cfg,header);
                            elseif isfield(cfg,'EXTRA') & isfield(cfg.EXTRA,'exclude')
                                cfg.exclude = cfg.EXTRA.exclude;
                            else
                                cfg.exclude = [];
                            end
                        end
                    end
                    
                    % Save settings
                    if skipprocessing == 0 & ~exist('flagsettings','var') & cfg.SettingsFileAuto == false
                        lab_save_settings(cfg,calc);
                        if isfield(cfg,'IS') & isfield(cfg.IS,'LF')
                            if isfield(cfg.IS.LF,'correctdigits') & cfg.IS.LF.correctdigits == true
                                lab_correct_digits(calc.Filelist);
                            end
                            if isfield(cfg.IS.LF,'coregmode') & cfg.IS.LF.coregmode > 1
                                lab_prepare_coreg(calc.Filelist,cfg.IS);
                            end
                        elseif isfield(cfg,'ISFFT') & isfield(cfg.ISFFT,'LF')
                            if isfield(cfg.ISFFT.LF,'correctdigits') & cfg.ISFFT.LF.correctdigits == true
                                lab_correct_digits(calc.Filelist);
                            end
                            if isfield(cfg.ISFFT.LF,'coregmode') & cfg.ISFFT.LF.coregmode > 1
                                lab_prepare_coreg(calc.Filelist,cfg.ISFFT);
                            end
                        end
                        flagsettings = 1; %#ok<NASGU>
                    end
                    
                    % Preprocessing
                    if skipprocessing == 0 & skippreprocessing == 0
                        if size(calc.Filelist,2) == 1 & cfg.MAIN.auto == 0
                            [data,header,cfg,skipprocessing,skiperror] = lab_preprocessing(data,header,cfg,calc,skipfilt);
                        else
                            try
                                [data,header,cfg,skipprocessing,skiperror] = lab_preprocessing(data,header,cfg,calc,skipfilt);
                            catch err
                                disp(getReport(err))
                                skipprocessing = 1;
                                skiperror = 1;
                            end
                        end
                    elseif skipprocessing == 0 & isfield(cfg,'STITCHALL') & ~isempty(cfg.STITCHALL)
                        [data,header,cfg,skipprocessing] = lab_stitchingall(data,header,cfg);
                    end
                    
                    % Postprocessing
                    if skipprocessing == 0 & skippostprocessing == 0
                        if size(calc.Filelist,2) == 1 & cfg.MAIN.auto == 0
                            cfg = lab_postprocessing(data,header,cfg);
                        else
                            try
                                cfg = lab_postprocessing(data,header,cfg);
                            catch err
                                disp(getReport(err))
                                skipprocessing = 1;
                                skiperror = 1;
                            end
                        end
                    end
                    
                    if skiperror == 0
                        calc.Filelist_done = union(calc.Filelist_done,cellstr(WorkingEEG_file));
                        if StoreAll == 1
                            calc.Filelist_doneall = union(calc.Filelist_doneall,calc.Filelist(1,filenr));
                        end
                        % Write Filelists
                        if isfield(cfg,'SEARCH') & isfield(cfg.SEARCH,'searchfolder') & ~isempty(cfg.SEARCH.searchfolder) & ...
                                cfg.SEARCH.searchfolder ~= 0 & isfield(cfg,'settings_path') & ~isempty(cfg.settings_path)
                            if ~isempty(calc.Filelist_doneall)
                                if ispc
                                    lab_write_xls(fullfile(cfg.settings_path,'ProcessedFiles.xls'),calc.Filelist_doneall');
                                else
                                    if exist(fullfile(cfg.settings_path,'ProcessedFiles.txt'),'file')
                                        delete(fullfile(cfg.settings_path,'ProcessedFiles.txt'));
                                    end
                                    lab_write_txt(fullfile(cfg.settings_path,'ProcessedFiles.txt'),calc.Filelist_doneall');
                                end
                            end
                            if ~isempty(calc.Filelist_done)
                                if ispc
                                    lab_write_xls(fullfile(cfg.settings_path,'ProcessedFilesAll.xls'),calc.Filelist_done');
                                else
                                    if exist(fullfile(cfg.settings_path,'ProcessedFilesAll.txt'),'file')
                                        delete(fullfile(cfg.settings_path,'ProcessedFilesAll.txt'));
                                    end
                                    lab_write_txt(fullfile(cfg.settings_path,'ProcessedFilesAll.txt'),calc.Filelist_done');
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end
if skipprocessing == 0 & cfg.SettingsFileAuto == false
    % Save settings (cfg) to settings.mat
    lab_save_settings(cfg,calc);
end
disp('Finished processing files')

% Turn of log-file
diary off;
