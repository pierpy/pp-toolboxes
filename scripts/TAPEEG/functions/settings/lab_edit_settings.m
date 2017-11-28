% written by F. Hatz 2012

function cfg = lab_edit_settings(settingsfile,mode)

cfg = [];
nostore = 0;

if ~exist('settingsfile','var') | isempty(settingsfile)
    if ~exist('mode','var')
        mode = questdlg('Settings-file','Settings file','Read','Edit','Create','Edit');
    end
    if strcmp(mode,'Edit')
        disp ('Ask for settings file')
        [Settings_file,Settings_filepath] = uigetfile('*.mat','Select settings file');
        if Settings_file ~= 0
            settingsfile = fullfile(Settings_filepath,Settings_file);
            Title = ['Edit ' settingsfile];
        else
            return
        end
        clearvars Settings_file
    elseif strcmp(mode,'Create') | strcmp(mode,'CreateCollectGraph') |  strcmp(mode,'CreateCollectSpectra') | ...
            strcmp(mode,'CreateCollectConnectivity') |  strcmp(mode,'CreateCollectMicrostates') | ...
            strcmp(mode,'CreateMRI') | strcmp(mode,'CreateMatrix')
        [Settings_file,Settings_filepath] = uiputfile('*.mat','Select directory and file to store settings');
        if Settings_filepath ~= 0
            settingsfile = fullfile(Settings_filepath,Settings_file);
            if exist(settingsfile,'file')
                delete(settingsfile)
            end
            Title = ['Create ' settingsfile];
        else
            return
        end
        clearvars Settings_file
    elseif strcmp(mode,'Read')
        disp ('Ask for settings file')
        [Settings_file,Settings_filepath]=uigetfile('*.mat','Select settings file');
        if Settings_file ~= 0
            settingsfile = fullfile(Settings_filepath,Settings_file);
            nostore = 1;
            Title = ['Read ' settingsfile];
        else
            return
        end
        clearvars Settings_file
    elseif ischar(mode) & exist(mode,'dir')
        Settings_filepath = mode;
        settingsfile = fullfile(Settings_filepath,'settingsauto.mat');
        Title = ['Create Autosettings (' mode ')'];
    else
        return
    end
end

if exist(settingsfile,'file')
    disp('Load settings.mat (from previous run)')
    load(settingsfile);
end

if strcmp(mode,'CreateCollectSpectra') | (isfield(cfg,'CollectFFT') & ~isempty(cfg.CollectFFT))
    [cfg,~,skipprocessing] = lab_set_collect_spectraldata(cfg,true);
    if skipprocessing == 1
        return
    end
elseif strcmp(mode,'CreateCollectGraph') | (isfield(cfg,'CollectGraph') & ~isempty(cfg.CollectGraph))
    [cfg,~,skipprocessing] = lab_set_collect_graphanalysis(cfg,true);
    if skipprocessing == 1
        return
    end
elseif strcmp(mode,'CreateCollectConnectivity') | (isfield(cfg,'CollectConnect') & ~isempty(cfg.CollectConnect))
    [cfg,~,skipprocessing] = lab_set_collect_connectivity(cfg,true);
    if skipprocessing == 1
        return
    end
elseif strcmp(mode,'CreateCollectMicrostates') | (isfield(cfg,'CollectMicro') & ~isempty(cfg.CollectMicro))
    [cfg,~,skipprocessing] = lab_set_collect_microstates(cfg,true);
    if skipprocessing == 1
        return
    end
elseif strcmp(mode,'CreateMRI') | (isfield(cfg,'MRI') & ~isempty(cfg.MRI))
    [cfg,skipprocessing] = lab_set_process_mri(cfg,true);
    if skipprocessing == 1
        return
    end
elseif strcmp(mode,'CreateMatrix') | (isfield(cfg,'MATRIX') & ~isempty(cfg.MATRIX))
    [cfg,skipprocessing] = lab_set_process_matrix(cfg,false,true,false,true);
    if skipprocessing == 1
        return
    end
else
    cfg.SEARCH.searchfolder = Settings_filepath;
    [cfg,~,skipprocessing] = lab_edit_cfg(cfg,[],0,0,0,Title);
    if skipprocessing == 1
        return
    end
end

% Store settings
if nostore == 0
    disp('Save new settings')
    [settings_file,settings_path] = lab_filename(settingsfile);
    cfg.settings_file = settings_file;
    cfg.settings_path = settings_path;
    lab_save_settings(cfg);
else
    disp('Settings not saved, changes are lost')
end

return