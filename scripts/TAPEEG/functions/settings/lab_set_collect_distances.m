function [cfg,Files,skipprocessing] = lab_set_collect_distances(cfg,nofolder)

skipprocessing = 0;
Files = {};

if ~exist('nofolder','var')
    nofolder = false;
end
if ~exist('cfg','var') | ~isfield(cfg,'CollectDistance') | ~isfield(cfg.CollectDistance,'searchfolder')
    cfg.CollectDistance.searchfolder = '';
    cfg.CollectDistance.includestring{1} = '';
    cfg.CollectDistance.excludestring{1} = '';
    cfg.CollectDistance.subjecttext = ' ';
    cfg.CollectDistance.subjectname = [];
    cfg.CollectDistance.outputfolder = 'Distances';
    cfg.CollectDistance.Threshold = 5;
end
if  nofolder == true
    cfg.CollectDistance.searchfolder = '';
    cfg.CollectDistance.Files = {};
elseif ~isfield(cfg.CollectDistance,'Files')
    cfg.CollectDistance.Files = {};
end

Formats = {};
Prompt = cell(0,2);

if nofolder == false
    Prompt{end+1,1} = 'Search folder';
    Formats(end+1,1).type = 'text';
    Formats(end,1).span = [1 2];
    
    Prompt(end+1,:) = {'','searchfolder'};
    Formats(end+1,1).type = 'edit';
    Formats(end,1).format = 'dir';
    Formats(end,1).size = 300;
    Formats(end,1).span = [1 2];
    Formats(end,1).callback = {@set_folder,'@ALL','@ALL'};
end

Prompt{end+1,1} = 'Include strings';
Formats(end+1,1).type = 'text';
Formats(end,1).span = [1 2];

Prompt(end+1,:) = {'','includestring'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'vector';
Formats(end,1).size = 300;
Formats(end,1).span = [1 2];

Formats(end+1,1).type = 'none';
Formats(end,1).span = [1 2];

Prompt{end+1,1} = 'Exclude strings';
Formats(end+1,1).type = 'text';
Formats(end,1).span = [1 2];

Prompt(end+1,:) = {'','excludestring'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'vector';
Formats(end,1).size = 300;
Formats(end,1).span = [1 2];

Prompt{end+1,1} = '';
Formats(end+1,1).type = 'text';
Formats(end,1).span = [1 2];

if nofolder == false
    Prompt(end+1,:) = {'Search Files','Files'};
    Formats(end+1,1).type = 'check';
    Formats(end,1).callback = {@search_files,'@ALL','@ALL','$subjecttext'};
end

Prompt(end+1,:) = {'Preload File-Info',''};
Formats(end+1,1).type = 'button';
Formats(end,1).style = 'pushbutton';
Formats(end,1).size = [110 25];
Formats(end,1).callback = {@load_fileinfo,'@ALL','@ALL','$subjecttext'};

Prompt(end+1,:) = {' ',''};
Formats(end+1,1).type = 'text';
Formats(end,1).span = [1 2];

Prompt(end+1,:) = {cfg.CollectDistance.subjecttext,'subjecttext'};
Formats(end+1,1).type = 'text';
Formats(end,1).span = [1 2];

Prompt(end+1,:) = {'Number of underscores in subject name','subjectname'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'integer';
Formats(end,1).limits = [-1 99];
Formats(end,1).size = 30;
Formats(end,1).span = [1 2];

Formats(end+1,1).type = 'none';
Formats(end,1).span = [1 2];

Prompt(end+1,:) = {'Output-folder', 'outputfolder'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'text';
Formats(end,1).size = 160;
Formats(end,1).span = [1 2];

Formats(end+1,1).type = 'none';
Formats(end,1).span = [1 2];

Prompt(end+1,:) = {'Threshold','Threshold'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'float';
Formats(end,1).limits = [0 inf];
Formats(end,1).size = 40;

[cfg.CollectDistance,Cancelled] = inputsdlg(Prompt,'Collect Distances',Formats,cfg.CollectDistance);
if Cancelled == 1
    cfg.CollectDistance = [];
    skipprocessing = 1;
else
    pause(0.2);
    Files = cfg.CollectDistance.Files;
    cfg.CollectDistance = rmfield(cfg.CollectDistance,'Files');
    if isfield(cfg.CollectDistance,'searchfolder') & exist(cfg.CollectDistance.searchfolder,'dir')
        save(fullfile(cfg.CollectDistance.searchfolder,'settings_distance.mat'),'cfg','-v7.3');
    end
end

end

function settings = set_folder(settings)
    if exist(fullfile(settings.searchfolder,'settings_distance.mat'),'file')
        searchfolder = settings.searchfolder;
        load(fullfile(settings.searchfolder,'settings_distance.mat'))
        if exist('cfg','var') & isfield(cfg,'CollectDistance') & ~isempty(cfg.CollectDistance)
            settings = cfg.CollectDistance;
            settings.searchfolder = searchfolder;
        end
    end
end

function settings = load_fileinfo(settings,Shandle)
    Files = settings.Files;
    if isempty(Files)
        [Filename,Filepath] = uigetfile('*.mat','Select *Distance.mat');
        if isnumeric(Filename)
            return
        end
        Files{1} = fullfile(Filepath,Filename);
    end
    Subjectname = {'',''};
    if ~isempty(Files);
        try %#ok<TRYNC>
            MAT = load(Files{1,1});
            Subjectname = lab_prepare_subjectname(Files{1,1});
            if isfield(MAT,'patient')
                Subjectname{1} = ['(' MAT.patient '^-^1)'];
                Subjectname{1} = regexprep(Subjectname{1},'_',' ');
                if isempty(settings.subjectname)
                    settings.subjectname = -1;
                end
            end
        end
    end
    if ~isempty(Subjectname{1})
        settings.subjecttext = [Subjectname{1} ' ' Subjectname{2}];
    else
        settings.subjecttext = Subjectname{2};
    end
    set(Shandle,'String',settings.subjecttext);
end

function settings = search_files(settings,Shandle)
    if isempty(settings.Files)
        cfg2.CollectDistance = settings;
        Files = lab_collect_distances_search(cfg2);
        if ~isempty(Files)
            settings.Files = Files;
        end
    else
        disp ('Select Files')
        selection = listdlg('PromptString','Files:','SelectionMode','multiple', ...
            'ListString',settings.Files,'InitialValue',1:length(settings.Files),'CancelString','None','ListSize',[450 400]);
        pause(0.2);
        if ~isempty(selection)
            settings.Files = settings.Files(1,selection);
        else
            settings.Files = {};
        end
    end
    if ~isempty(settings.Files)
        settings = load_fileinfo(settings,Shandle);
    end
end