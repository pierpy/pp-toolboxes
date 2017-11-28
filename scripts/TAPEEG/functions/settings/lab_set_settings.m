function [cfg,skipprocessing] = lab_set_settings(cfg)

global Main_Path

skipprocessing = 0;

if ~exist('cfg','var') | ~isfield(cfg,'settings_path')
    cfg.settings_path = pwd;
end
disp(['Search/Open settings-files (' cfg.settings_path ')']);

% Load settings.mat from previous run
Strings = {'settings*.mat','|settings_graph.mat','|settings_connect.mat','|settings_microstates.mat', ...
    '|settings_spectras.mat','|settings_connectivity.mat'};

if exist(fullfile(Main_Path,'settings'),'dir')
    [List] = lab_search(fullfile(Main_Path,'settings'),Strings,1,1);
    Dates = zeros(1,length(List));
else
    List = {};
    Dates = [];
end
if exist(cfg.settings_path,'dir')
    [List2,Param] = lab_search(cfg.settings_path,Strings,1,1);
    List = [List(:)' List2(:)'];
    if exist(fullfile(cfg.settings_path,'settings'),'dir')
        [List2,Param2] = lab_search(fullfile(cfg.settings_path,'settings'),Strings,1);
        List = [List(:)' List2(:)'];
        Param = [Param(:)' Param2(:)'];
        clearvars List2 Param2
    end
    if ~isempty(List)
        Param = struct2cell(Param);
        Param = permute(Param,[3 1 2]);
        Dates2 = cell2mat(Param(:,5));
        clearvars Param
        Dates = [Dates Dates2(:)'];
    end
end
clearvars Strings
if isfield(cfg,'MAIN') & isfield(cfg.MAIN,'auto') & cfg.MAIN.auto == 1
    if isempty(List)
        skipprocessing = 1;
        return
    end
    if isfield(cfg,'settings_fileM') & exist(cfg.settings_fileM,'file')
        tmp = find(strcmp(List,cfg.settings_fileM),1);
    else
        tmp = find(Dates == max(Dates),1,'first');
    end
    if ~isempty(tmp)
        settings_file = List{tmp};
    else
        settings_file = List{1};
    end
    clearvars tmp
else
    List = [cellstr('New settings') List];
    Dates = [0 Dates(:)'];
    if isfield(cfg,'settings_fileM') & exist(cfg.settings_fileM,'file')
        tmp = find(strcmp(List,cfg.settings_fileM),1);
    else
        tmp = find(Dates == max(Dates),1,'first');
    end
    if ~isempty(tmp)
        cfg.selection = List{tmp};
    elseif size(List,2) > 1
        cfg.selection = List{end};
    else
        cfg.selection = List{1};
    end
    clearvars tmp
    cfg.editsettings = true;
    
    Prompt = cell(0,2);
    Formats = [];

    Prompt(end+1,:) = {'Settings-File','selection'};
    Formats(end+1,1).type = 'List';
    Formats(end,1).style = 'popupmenu';
    Formats(end,1).format = 'input';
    Formats(end,1).size = 550;
    Formats(end,1).items = List;
    
    Prompt(end+1,:) = {'Edit','editsettings'};
    Formats(end+1,1).type = 'check';
    
    Prompt(end+1,:) = {'Search',''};
    Formats(end+1,1).type = 'button';
    Formats(end,1).style = 'pushbutton';
    Formats(end,1).size = [100 25];
    Formats(end,1).callback = {@search_settings,'@ALL','@ALL'};

    [cfg,Cancelled] = inputsdlg(Prompt,'Load settings',Formats,cfg);
    if Cancelled == 1
        cfg = rmfield(cfg,'selection');
        cfg.editsettings = true;
        return
    else
        pause(0.2);
        settings_file = cfg.selection;
        cfg = rmfield(cfg,'selection');
    end
end

if exist(settings_file,'file')
    disp(['Load ' lab_filename(settings_file) ' (from previous run)'])
    cfgorig = cfg;
    load(settings_file);
    if isfield(cfgorig,'MAIN')
        cfg.MAIN.auto = cfgorig.MAIN.auto;
        if isfield(cfgorig.MAIN,'log_path')
            cfg.MAIN.log_path = cfgorig.MAIN.log_path;
        end
    end
    if isfield(cfgorig,'SEARCH') & isfield(cfgorig.SEARCH,'searchfolder')
        cfg.SEARCH.searchfolder = cfgorig.SEARCH.searchfolder;
    end
    cfg.settings_path = cfgorig.settings_path;
    if isfield(cfgorig,'settings_file')
        cfg.settings_file = cfgorig.settings_file;
    else
        cfg.settings_file = lab_filename(settings_file);
    end
    if isfield(cfgorig,'editsettings')
        cfg.editsettings = cfgorig.editsettings;
    end
else
    cfg.editsettings = true;
    SettingsName = inputdlg('Settings-Name');
    cfg.settings_file = [SettingsName{1,1} '_settings.mat'];
end

end

function settings = search_settings(settings)
   [Settings_file,Settings_filepath] = uigetfile('*.mat','Select settings-file');
   if ~isnumeric(Settings_file)
       settings.selection = fullfile(Settings_filepath,Settings_file);
   end
end