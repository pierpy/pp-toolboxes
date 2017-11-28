function [cfg,skipprocessing] = lab_set_process_matrix(cfg,nofolder,setsearch,nograph,donumber,dowrite,header)

global MatrixAll

skipprocessing = 0;
if ~exist('header','var')
    header = [];
end
if ~exist('dowrite','var')
    dowrite = true;
end
if ~exist('donumber','var')
    donumber = false;
end
if ~exist('nograph','var')
    nograph = false;
end
if ~exist('setsearch','var')
    setsearch = false;
end
if ~exist('nofolder','var')
    nofolder = false;
end

if ~exist('cfg','var') | ~isfield(cfg,'MATRIX') | ~isfield(cfg.MATRIX,'outputfolder')
    cfg.MATRIX.nummatrices = [];
    cfg.MATRIX.exclude = [];
    cfg.MATRIX.donormalize = 'off';
    cfg.MATRIX.dosymmetrical = false;
    cfg.MATRIX.dodirectional = false;
    cfg.MATRIX.minzero = true;
    cfg.MATRIX.rankmatrix = false;
    cfg.MATRIX.rankorder = 5;
    cfg.MATRIX.grandaverage = false;
    cfg.MATRIX.Mappings = [];
    cfg.MATRIX.MappingsMode = 'Average';
    cfg.MATRIX.MappingsWrite = false;
    cfg.MATRIX.binary = false;
    cfg.MATRIX.binarythreshold = 0.5;
    cfg.MATRIX.binarymode = 'fixed';
    cfg.MATRIX.binarymethod = 'Threshold';
    if nofolder == false
        cfg.MATRIX.outputfolder = 'MatrixResults';
    else
        cfg.MATRIX.outputfolder = '';
    end
    cfg.MATRIX.SEARCH.searchstring{1,1} = 'matrix.txt';
end

if isfield(cfg,'Output_file')
    Subjectname = lab_prepare_subjectname(fullfile(cfg.Output_filepath,cfg.Output_file));
    if ~isempty(Subjectname{1})
        cfg.MATRIX.subjecttext = [Subjectname{1} ' ' Subjectname{2}];
    else
        cfg.MATRIX.subjecttext = Subjectname{2};
    end
    if ~isfield(cfg.MATRIX,'subjectname')
        if isfield(cfg,'subjectname') & cfg.subjectname >= 0
            cfg.MATRIX.subjectname = cfg.subjectname;
        else
            cfg.MATRIX.subjectname = 0;
        end
    end
else
    cfg.MATRIX.subjecttext = '';
end

if setsearch == true & isfield(cfg,'SEARCH')
    cfg.MATRIX.SEARCH = cfg.SEARCH;
end

Prompt = {};
Formats = {};

if setsearch == true
    Prompt(end+1,:) = {'Search Files','SEARCH'};
    Formats(end+1,1).type = 'check';
    Formats(end,1).callback = {@lab_set_searchstrings,'@ALL','@ALL',0,0,{'strings','Matrix (.txt)','Matrix (.mat)'}};
    Formats(end,1).span = [1 4];
    
    Formats(end+1,1).type = 'none';
    Formats(end,1).span = [1 4];
end

if setsearch == true | isfield(cfg,'Output_file')
    Prompt(end+1,:) = {cfg.MATRIX.subjecttext,'subjecttext'};
    Formats(end+1,1).type = 'text';
    Formats(end,1).span = [1 4];
    
    Prompt(end+1,:) = {'Number of underscores in subject name','subjectname'};
    Formats(end+1,1).type = 'edit';
    Formats(end,1).format = 'integer';
    Formats(end,1).limits = [-1 99];
    Formats(end,1).size = 30;
    Formats(end,1).span = [1 3];
    
    Prompt(end+1,:) = {'Preload File-Info',''};
    Formats(end+1,1).type = 'button';
    Formats(end,1).style = 'pushbutton';
    Formats(end,1).size = [110 25];
    Formats(end,1).callback = {@load_fileinfo,'@ALL','@ALL','$subjecttext'};
end

if nofolder == false
    Prompt(end+1,:) = {'Result folder', 'outputfolder'};
    Formats(end+1,1).type = 'edit';
    Formats(end,1).format = 'text';
    Formats(end,1).size = 100;
    Formats(end,1).span = [1 4];
end

Formats(end+1,1).type = 'none';
Formats(end,1).span = [1 4];

if donumber == true
    Prompt(end+1,:) = {'Number of matrices per subject (empty = all)','nummatrices'};
    Formats(end+1,1).type = 'edit';
    Formats(end,1).format = 'integer';
    Formats(end,1).limits = [0 9999999];
    Formats(end,1).size = 40;
    Formats(end,1).span = [1 3];
end

Prompt(end+1,:) = {'Exclude channels','exclude'};
Formats(end+1,1).type = 'check';
Formats(end,1).callback = {@lab_get_exclude,'exclude','exclude'};
Formats(end,1).span = [1 4];

Prompt(end+1,:) = {'Normalize matrix','donormalize'};
Formats(end+1,1).type = 'list';
Formats(end,1).style = 'popupmenu';
Formats(end,1).format = 'input';
Formats(end,1).items = {'off','single','over all'};
Formats(end,1).span = [1 4];

Prompt(end+1,:) = {'Transform to symmetrical matrix','dosymmetrical'};
Formats(end+1,1).type = 'check';
Formats(end,1).span = [1 4];

Prompt(end+1,:) = {'Transform to directional matrix (eg dPLI)','dodirectional'};
Formats(end+1,1).type = 'check';
Formats(end,1).span = [1 4];

Prompt(end+1,:) = {'Set values < 0 to zero','minzero'};
Formats(end+1,1).type = 'check';
Formats(end,1).span = [1 4];

Prompt(end+1,:) = {'Correct distance','DIST'};
Formats(end+1,1).type = 'check';
Formats(end,1).callback = {@lab_get_DIST,'DIST','DIST'};
Formats(end,1).span = [1 4];

Prompt(end+1,:) = {'Rank matrix','rankmatrix'};
Formats(end+1,1).type = 'check';
Formats(end,1).span = [1 2];

Prompt(end+1,:) = {'Order','rankorder'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'float';
Formats(end,1).limits = [0 99];
Formats(end,1).size = 40;
Formats(end,1).span = [1 2];

Prompt(end+1,:) = {'           ',''};
Formats(end+1,1).type = 'text';
Formats(end,1).span = [1 4];

if donumber == true
    Prompt(end+1,:) = {'Average','average'};
    Formats(end+1,1).type = 'list';
    Formats(end,1).style = 'popupmenu';
    Formats(end,1).format = 'input';
    Formats(end,1).items = {'off','subject','folder','all'};
    Formats(end,1).span = [1 4];
else
    Prompt(end+1,:) = {'Average','average'};
    Formats(end+1,1).type = 'list';
    Formats(end,1).style = 'popupmenu';
    Formats(end,1).format = 'input';
    Formats(end,1).items = {'off','folder','all'};
    Formats(end,1).span = [1 4];
end

Formats(end+1,1).type = 'none';
Formats(end,1).span = [1 4];

Formats(end+1,1).type = 'none';

Prompt(end+1,:) = {'Mappings','Mappings'};
Formats(end+1,1).type = 'check';
Formats(end,1).size = [-1 -1];
Formats(end,1).callback = {@lab_load_mappings,'Mappings','Mappings',cfg};

Prompt(end+1,:) = {'','MappingsMode'};
Formats(end+1,1).type = 'list';
Formats(end,1).style = 'popupmenu';
Formats(end,1).format = 'input';
Formats(end,1).items = {'Average','Degree','Betweenness'};

if dowrite == true
    Prompt(end+1,:) = {'xls-write','MappingsWrite'};
    Formats(end+1,1).type = 'check';
else
    Formats(end+1,1).type = 'none';
end

Formats(end+1,1).type = 'none';

Prompt(end+1,:) = {'Make binary','binary'};
Formats(end+1,1).type = 'check';

Prompt(end+1,:) = {'Threshold','binarythreshold'};
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

if nograph == false
    Formats(end+1,1).type = 'none';
    
    Prompt(end+1,:) = {'Plot matrices','PLOT'};
    Formats(end+1,1).type = 'check';
    Formats(end,1).callback = {@lab_set_plot_matrices,'PLOT','PLOT',[],'Mappings'};
    Formats(end,1).span = [1 3];
    
    Formats(end+1,1).type = 'none';
    
    Prompt(end+1,:) = {'Graph analysis','GRAPH'};
    Formats(end+1,1).type = 'check';
    Formats(end,1).callback = {@lab_set_graphanalysis,'@ALL','@ALL',[],cfg,header};
    Formats(end,1).span = [1 3];
end

[cfg.MATRIX,Cancelled] = inputsdlg(Prompt,'Process matrix',Formats,cfg.MATRIX);
if Cancelled == 1
    skipprocessing = 1;
    cfg.MATRIX = [];
    return
else
    MatrixAll = [];
    pause(0.2);
end

if setsearch == true & ~isempty(cfg.MATRIX.SEARCH)
    cfg.SEARCH = cfg.MATRIX.SEARCH;
    cfg.MATRIX = rmfield(cfg.MATRIX,'SEARCH');
end

    function settings = load_fileinfo(settings,Shandle)
        [Filename,Filepath] = uigetfile('*.mat','Select *Graph.mat');
        if isnumeric(Filename)
            return
        end
        Filename = fullfile(Filepath,Filename);
        Subjectname = {'',''};
        if ~isnumeric(Filename) & ~isempty(Filename)
            Subjectname = lab_prepare_subjectname(Filename);
            try %#ok<TRYNC>
                [~,~,cfg] = lab_read_matrix(Filename,cfg,false,true);
                if isfield(cfg,'patient')
                    Subjectname{1} = ['(' cfg.patient '^-^1)'];
                    Subjectname{1} = regexprep(Subjectname{1},'_',' ');
                    if isempty(settings.subjectname)
                        settings.subjectname = -1;
                    end
                end
                clearvars cfg2
            end
        end
        if ~isempty(Subjectname{1})
            settings.subjecttext = [Subjectname{1} ' ' Subjectname{2}];
        else
            settings.subjecttext = Subjectname{2};
        end
        set(Shandle,'String',settings.subjecttext);
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
    
    function exclude = lab_get_exclude(exclude)
        settings.exclude = exclude;
        Prompt2 = {'Channels','exclude'};
        Formats2.type = 'edit';
        Formats2.format = 'vector';
        Formats2.size = 150;
        Formats2.limits = [-inf inf];
        [settings,Cancelled2] = inputsdlg(Prompt2,'Exclude channels',Formats2,settings);
        if Cancelled2 == 1
            exclude = [];
        else
            exclude = settings.exclude;
        end
    end
end