function [cfg,Files,skipprocessing] = lab_set_collect_connectivity(cfg,nofolder)

skipprocessing = 0;
Files = [];
Patient = ' ';
Filename = ' ';
Header = [];

if ~exist('nofolder','var')
    nofolder = false;
end
if ~exist('cfg','var') | ~isfield(cfg,'CollectConnect') | ~isfield(cfg.CollectConnect,'searchfolder')
    cfg.CollectConnect.searchfolder = '';
    cfg.CollectConnect.includestring{1} = '';
    cfg.CollectConnect.excludestring{1} = '';
    cfg.CollectConnect.Files = [];
    cfg.CollectConnect.subjecttext = ' ';
    cfg.CollectConnect.subjectname = [];
    cfg.CollectConnect.SelectVars = false;
    cfg.CollectConnect.outputfolder = 'ConnectivityAnalysis';
    cfg.CollectConnect.randphase = false;
    cfg.CollectConnect.valuerandphase = [];
    cfg.CollectConnect.value2randphase = ' ';
    cfg.CollectConnect.methodrandphase = 'P-Value';
    cfg.CollectConnect.binary = false;
    cfg.CollectConnect.binarythreshold = 0.5;
    cfg.CollectConnect.binarymode = 'fixed';
    cfg.CollectConnect.doaverage = false;
    cfg.CollectConnect.nummatrices = [];
    cfg.CollectConnect.Matrix = [];
    cfg.CollectConnect.WriteMatrices = true;
    cfg.CollectConnect.WriteDegrees = true;
    cfg.CollectConnect.WriteConnections = false;
    cfg.CollectConnect.PLOT = [];
    cfg.CollectConnect.GRAPH = [];
    cfg.CollectConnect.Kmeans = [];
end
if  nofolder == true
    cfg.CollectConnect.searchfolder = '';
    cfg.CollectConnect.Files = [];
elseif ~isfield(cfg.CollectConnect,'Files')
    cfg.CollectConnect.Files = [];
end

Formats = {};
Prompt = cell(0,2);

if nofolder == false
    Prompt{end+1,1} = 'Search folder';
    Formats(end+1,1).type = 'text';
    Prompt(end+1,:) = {'','searchfolder'};
    Formats(end+1,1).type = 'edit';
    Formats(end,1).format = 'dir';
    Formats(end,1).size = 300;
    Formats(end,1).span = [1 5];
    Formats(end,1).callback = {@set_folder,'@ALL','@ALL'};
end

Prompt{end+1,1} = 'Include strings';
Formats(end+1,1).type = 'text';
Prompt(end+1,:) = {'','includestring'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'vector';
Formats(end,1).size = 300;
Formats(end,1).span = [1 5];

Prompt{end+1,1} = '';
Formats(end+1,1).type = 'text';
Formats(end,1).span = [1 5];

Prompt{end+1,1} = 'Exclude strings';
Formats(end+1,1).type = 'text';
Prompt(end+1,:) = {'','excludestring'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'vector';
Formats(end,1).size = 300;
Formats(end,1).span = [1 5];

Prompt{end+1,1} = '';
Formats(end+1,1).type = 'text';
Formats(end,1).span = [1 5];

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
Formats(end,1).span = [1 4];

Formats(end+1,1).type = 'none';
Formats(end,1).span = [1 5];

Prompt(end+1,:) = {cfg.CollectConnect.subjecttext,'subjecttext'};
Formats(end+1,1).type = 'text';
Formats(end,1).span = [1 5];

Prompt(end+1,:) = {'Number of underscores in subject name','subjectname'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'integer';
Formats(end,1).limits = [-99 99];
Formats(end,1).size = 30;
Formats(end,1).callback = {@set_patient,'@ALL','@ALL','$patienttext'};
Formats(end,1).span = [1 3];

Prompt(end+1,:) = {Patient,'patienttext'};
Formats(end+1,1).type = 'text';

Formats(end+1,1).type = 'none';
Formats(end,1).span = [1 5];

Prompt(end+1,:) = {'Output-folder', 'outputfolder'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'text';
Formats(end,1).size = 160;
Formats(end,1).span = [1 3];

Prompt(end+1,:) = {'Select connectivity measures','SelectVars'};
Formats(end+1,1).type = 'check';
Formats(end,1).span = [1 2];

Formats(end+1,1).type = 'none';
Formats(end,1).span = [1 5];

Prompt(end+1,:) = {'Correct by phase randomization','randphase'};
Formats(end+1,1).type = 'check';
Formats(end,1).span = [1 2];

Prompt(end+1,:) = {'','valuerandphase'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'float';
Formats(end,1).limits = [0 100];
Formats(end,1).size = 30;

Prompt(end+1,:) = {'','value2randphase'};
Formats(end+1,1).type = 'list';
Formats(end,1).style = 'popupmenu';
Formats(end,1).format = 'input';
Formats(end,1).items = {' ','std','percent'};
Formats(end,1).callback = {@set_randphase,'@ALL','@ALL'};

Prompt(end+1,:) = {'Method','methodrandphase'};
Formats(end+1,1).type = 'list';
Formats(end,1).style = 'popupmenu';
Formats(end,1).format = 'input';
Formats(end,1).items = {'P-Value','Diff','Threshold'};
Formats(end,1).callback = {@set_method,'@ALL','@ALL'};

Formats(end+1,1).type = 'none';
Formats(end,1).span = [1 5];

Prompt(end+1,:) = {'Make binary','binary'};
Formats(end+1,1).type = 'check';

Prompt(end+1,:) = {'','binarythreshold'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'float';
Formats(end,1).limits = [-inf inf];
Formats(end,1).size = 40;

Prompt(end+1,:) = {'','binarymode'};
Formats(end+1,1).type = 'list';
Formats(end,1).style = 'popupmenu';
Formats(end,1).format = 'input';
Formats(end,1).items = {'fixed','std','percent'};
Formats(end,1).callback = {@set_binarymode,'@ALL','@ALL'};
Formats(end,1).span = [1 3];

Formats(end+1,1).type = 'none';
Formats(end,1).span = [1 5];

Prompt(end+1,:) = {'Number of matrices per subject (empty = all)','nummatrices'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'integer';
Formats(end,1).limits = [1 inf];
Formats(end,1).size = 40;
Formats(end,1).span = [1 5];

Prompt(end+1,:) = {'Average connectivity matrix per subject','doaverage'};
Formats(end+1,1).type = 'check';
Formats(end,1).span = [1 5];

Formats(end+1,1).type = 'none';
Formats(end,1).span = [1 5];

Prompt(end+1,:) = {'Process matrices','MATRIX'};
Formats(end+1,1).type = 'check';
Formats(end,1).callback = {@lab_set_process_matrix,'@ALL','@ALL',true,false,true,false,false};
Formats(end,1).span = [1 5];

Formats(end+1,1).type = 'none';
Formats(end,1).span = [1 5];

Prompt(end+1,:) = {'Write matrices','WriteMatrices'};
Formats(end+1,1).type = 'check';

Prompt(end+1,:) = {'Write degrees','WriteDegrees'};
Formats(end+1,1).type = 'check';
Formats(end,1).span = [1 2];

Prompt(end+1,:) = {'Write connections','WriteConnections'};
Formats(end+1,1).type = 'check';

Formats(end+1,1).type = 'none';
Formats(end,1).span = [1 5];

Prompt(end+1,:) = {'Plot matrices','PLOT'};
Formats(end+1,1).type = 'check';
Formats(end,1).callback = {@lab_set_plot_matrices,'PLOT','PLOT'};

Prompt(end+1,:) = {'Graph analysis','GRAPH'};
Formats(end+1,1).type = 'check';
Formats(end,1).callback = {@set_graphanalysis,'@ALL','@ALL'};
Formats(end,1).span = [1 2];

Prompt(end+1,:) = {'Kmeans clustering','Kmeans'};
Formats(end+1,1).type = 'check';
Formats(end,1).callback = {@lab_set_calculate_kmeans,'Kmeans','Kmeans',0};

[cfg.CollectConnect,Cancelled] = inputsdlg(Prompt,'Collect Connectivity',Formats,cfg.CollectConnect);
if Cancelled == 1
    cfg.CollectConnect = [];
    skipprocessing = 1;
else
    pause(0.2);
    Files = cfg.CollectConnect.Files;
    cfg.CollectConnect = rmfield(cfg.CollectConnect,'Files');
    if isfield(cfg.CollectConnect,'searchfolder') & exist(cfg.CollectConnect.searchfolder,'dir')
        save(fullfile(cfg.CollectConnect.searchfolder,'settings_connectivity.mat'),'cfg','-v7.3');
    end
end

    function settings = load_fileinfo(settings,Shandle,Phandle)
        Filestmp = settings.Files;
        if isempty(Filestmp) | ~isfield(Filestmp,'list') | isempty(Filestmp(1).list)
            [Filename,Filepath] = uigetfile('*.mat','Select Conn_F...mat');
            if isnumeric(Filename)
                return
            end
            Filestmp.list{1} = fullfile(Filepath,Filename);
        end
        Subjectname = {'',''};
        if ~isempty(Filestmp) & isfield(Filestmp,'list') & ~isempty(Filestmp(1).list)
            Filename = Filestmp(1).list{1};
            try %#ok<TRYNC>
                MAT = load(Filename);
                Subjectname = lab_prepare_subjectname(Filename);
                if isfield(MAT,'patient')
                    Patient = MAT.patient;
                    set(Phandle,'String',['= ' regexprep(Patient,'_',' ')]);
                    if isempty(settings.subjectname)
                        settings.subjectname = [];
                    end
                end
                if isfield(MAT,'result') & isfield(MAT.result,'locs')
                    Header.locs = MAT.result.locs;
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
            cfgtmp.CollectConnect = settings;
            Filestmp = lab_collect_connectivity_search(cfgtmp);
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
                    [Ftmp,~,Itmp] = intersect(FilesAll,settings.Files(i).list);
                    if ~isempty(Ftmp)
                        settings.Files(i).list = Ftmp(:)';
                        if isfield(settings.Files,'listrand')
                            settings.Files(i).listrand = settings.Files(i).listrand(1,Itmp);
                        end
                        clearvars Ftmp Itmp
                        flags(i) = 1;
                    end
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
    
    function settings = set_randphase(settings)
        if strcmp(settings.value2randphase,'std')
            settings.valuerandphase = 1.28;
        elseif strcmp(settings.value2randphase,'percent')
            settings.valuerandphase = 80;
        end
        if strcmp(settings.value2randphase,' ')
            settings.valuerandphase = [];
            settings.methodrandphase = 'P-Value';
        elseif strcmp(settings.methodrandphase,'P-Value')
            settings.methodrandphase = 'Threshold';
        end
    end
    
    function settings = set_method(settings)
        if strcmp(settings.methodrandphase,'P-Value')
            settings.valuerandphase = [];
            settings.value2randphase = ' ';
        else
            settings.valuerandphase = 80;
            settings.value2randphase = 'percent';
        end
    end
    
    function settings = set_binarymode(settings)
        if strcmp(settings.binarymode,'std')
            settings.binarythreshold = 1.28;
        elseif strcmp(settings.binarymode,'percent')
            settings.binarythreshold = 80;
        else
            settings.binarythreshold = 0.5;
        end
    end
    
    function settings = set_graphanalysis(settings)
        settings = lab_set_graphanalysis(settings,[],cfg,Header,true);
    end

end

function settings = set_folder(settings)
    if exist(fullfile(settings.searchfolder,'settings_connectivity.mat'),'file')
        searchfolder = settings.searchfolder;
        load(fullfile(settings.searchfolder,'settings_connectivity.mat'))
        if exist('cfg','var') & isfield(cfg,'CollectConnect') & ~isempty(cfg.CollectConnect)
            settings = cfg.CollectConnect;
            settings.searchfolder = searchfolder;
        end
    end
end

