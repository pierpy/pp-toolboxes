function [cfg,Files,skipprocessing] = lab_set_collect_graphanalysis(cfg,nofolder)

skipprocessing = 0;
Files = [];
Patient = ' ';
Filename = ' ';

if ~exist('nofolder','var')
    nofolder = false;
end
if ~exist('cfg','var') | ~isfield(cfg,'CollectGraph') | ~isfield(cfg.CollectGraph,'searchfolder')
    cfg.CollectGraph.searchfolder = '';
    cfg.CollectGraph.includestring{1} = '';
    cfg.CollectGraph.excludestring{1} = '';
    cfg.CollectGraph.subjecttext = ' ';
    cfg.CollectGraph.subjectname = 0;
    cfg.CollectGraph.outputfolder = 'GraphAnalysis';
    cfg.CollectGraph.mode = 1;
    cfg.CollectGraph.nummatrices = [];
end
if  nofolder == true
    cfg.CollectGraph.searchfolder = '';
    cfg.CollectGraph.Files = [];
elseif ~isfield(cfg.CollectGraph,'Files')
    cfg.CollectGraph.Files = [];
end
if isempty(cfg.CollectGraph.mode)
    cfg.CollectGraph.mode = 1;
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
    Formats(end,1).size = 350;
    Formats(end,1).span = [1 3];
    Formats(end,1).callback = {@set_folder,'@ALL','@ALL'};
end

Prompt{end+1,1} = 'Include strings';
Formats(end+1,1).type = 'text';
Formats(end,1).span = [1 3];

Prompt(end+1,:) = {'','includestring'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'vector';
Formats(end,1).size = 350;
Formats(end,1).span = [1 3];

Formats(end+1,1).type = 'none';
Formats(end,1).span = [1 3];

Prompt{end+1,1} = 'Exclude strings';
Formats(end+1,1).type = 'text';
Formats(end,1).span = [1 3];

Prompt(end+1,:) = {'','excludestring'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'vector';
Formats(end,1).size = 350;
Formats(end,1).span = [1 3];

Formats(end+1,1).type = 'none';
Formats(end,1).span = [1 3];

if nofolder == false
    Prompt(end+1,:) = {'Search Files','Files'};
    Formats(end+1,1).type = 'check';
    Formats(end,1).callback = {@search_files,'@ALL','@ALL','$subjecttext','$patienttext'};
end

Prompt(end+1,:) = {'Preload File-Info',''};
Formats(end+1,1).type = 'button';
Formats(end,1).style = 'pushbutton';
Formats(end,1).size = [110 25];
Formats(end,1).callback = {@load_fileinfo,'@ALL','@ALL','$subjecttext','$patienttext'};

Formats(end+1,1).type = 'none';
Formats(end,1).span = [1 3];

Prompt(end+1,:) = {cfg.CollectGraph.subjecttext,'subjecttext'};
Formats(end+1,1).type = 'text';
Formats(end,1).span = [1 3];

Prompt(end+1,:) = {'Number of underscores in subject name','subjectname'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'integer';
Formats(end,1).limits = [-99 99];
Formats(end,1).size = 30;
Formats(end,1).callback = {@set_patient,'@ALL','@ALL','$patienttext'};
Formats(end,1).span = [1 2];

Prompt(end+1,:) = {Patient,'patienttext'};
Formats(end+1,1).type = 'text';

Formats(end+1,1).type = 'none';
Formats(end,1).span = [1 3];

Prompt(end+1,:) = {'Output-folder','outputfolder'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'text';
Formats(end,1).size = 160;
Formats(end,1).span = [1 3];

Formats(end+1,1).type = 'none';
Formats(end,1).span = [1 3];

Prompt(end+1,:) = {'','mode'};
Formats(end+1,1).type = 'list';
Formats(end,1).style = 'popupmenu';
Formats(end,1).items = {'Results of single matrices';'Average results of single matrices';'Results of average matrix'};
Formats(end,1).callback = {@set_mode,'@ALL','@ALL'};
Formats(end,1).span = [1 2];

Prompt(end+1,:) = {'Number of matrixes (empty = all)','nummatrices'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'integer';
Formats(end,1).limits = [1 inf];
Formats(end,1).size = 40;

[cfg.CollectGraph,Cancelled] = inputsdlg(Prompt,'Collect Graphanalysis',Formats,cfg.CollectGraph);
if Cancelled == 1
    cfg.CollectGraph = [];
    skipprocessing = 1;
else
    pause(0.2);
    Files = cfg.CollectGraph.Files;
    cfg.CollectGraph = rmfield(cfg.CollectGraph,'Files');
    if isfield(cfg.CollectGraph,'searchfolder') & exist(cfg.CollectGraph.searchfolder,'dir')
        save(fullfile(cfg.CollectGraph.searchfolder,'settings_graph.mat'),'cfg','-v7.3');
    end
end

    
    
    function settings = load_fileinfo(settings,Shandle,Phandle)
        Filestmp = settings.Files;
        if isempty(Filestmp) | ~isfield(Filestmp,'list') | isempty(Filestmp(1).list)
            [Filename,Filepath] = uigetfile('*.mat','Select *Graph.mat');
            if isnumeric(Filename)
                Filename = ' ';
                return
            end
            Filestmp.list{1} = fullfile(Filepath,Filename);
        end
        Subjectname = {'',''};
        if ~isempty(Filestmp) & isfield(Filestmp,'list') & ~isempty(Filestmp(1).list);
            Filename = Filestmp(1).list{1,1};
            Subjectname = lab_prepare_subjectname(Filename);
            try %#ok<TRYNC>
                MAT = load(Filename);
                if isfield(MAT,'Result') && isfield(MAT.Result,'patient')
                    Patient = MAT.Result.patient;
                    set(Phandle,'String',['= ' regexprep(Patient,'_',' ')]);
                    if isempty(settings.subjectname)
                        settings.subjectname = [];
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
    
    function settings = search_files(settings,Shandle,Phandle)
        if isempty(settings.Files)
            cfgtmp.CollectGraph = settings;
            Filestmp = lab_collect_graphanalysis_search(cfgtmp);
            if ~isempty(Filestmp)
                settings.Files = Filestmp;
            end
        else
            FilesAll = {};
            for i = 1:length(settings.Files)
                FilesAll = cat(2,FilesAll,settings.Files(i).list);
            end
            if isempty(FilesAll)
                settings.Files = [];
                return
            end
            disp ('Select Files')
            selection = listdlg('PromptString','Files:','SelectionMode','multiple', ...
                'ListString',FilesAll,'InitialValue',1:length(FilesAll),'CancelString','None','ListSize',[450 400]);
            pause(0.2);
            if ~isempty(selection)
                FilesAll = FilesAll(1,selection);
                flags = zeros(1,length(settings.Files));
                for i = 1:length(settings.Files)
                    Ftmp = intersect(FilesAll,settings.Files(i).list);
                    if ~isempty(Ftmp)
                        settings.Files(i).list = Ftmp(:)';
                        flags(i) = 1;
                    end
                    clearvars Ftmp
                end
                if max(flags) == 1
                    settings.Files = settings.Files(1,flags==1);
                else
                    settings.Files = [];
                end
                clearvars selection
            else
                settings.Files = [];
            end
        end
        if ~isempty(settings.Files) & isfield(settings.Files(1),'list') & ~isempty(settings.Files(1).list)
            settings = load_fileinfo(settings,Shandle,Phandle);
        end
    end
    
    function settings = set_patient(settings,Phandle)
        if ~isempty(settings.subjectname)
            set(Phandle,'String',['= ' regexprep(lab_subjectname(Filename,settings),'_',' ')]);
        else
            set(Phandle,'String',['= ' regexprep(Patient,'_',' ')]);
        end
    end
    
    function settings = set_mode(settings)
        if settings.mode == 3
            settings.nummatrices = [];
        end
    end

end

function settings = set_folder(settings)
    if exist(fullfile(settings.searchfolder,'settings_graph.mat'),'file')
        searchfolder = settings.searchfolder;
        load(fullfile(settings.searchfolder,'settings_graph.mat'))
        if exist('cfg','var') & isfield(cfg,'CollectGraph') & ~isempty(cfg.CollectGraph)
            settings = cfg.CollectGraph;
            settings.searchfolder = searchfolder;
            if isempty(settings.mode)
                settings.mode = 1;
            end
        end
    end
end