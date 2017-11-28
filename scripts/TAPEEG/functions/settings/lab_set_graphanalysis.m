function [settings,skipprocessing] = lab_set_graphanalysis(settings,nummatrices,cfg,header,nonumbers)

skipprocessing = 0;

if ~exist('nonumbers','var')
    nonumbers = false;
end
if ~exist('header','var')
    header = [];
end
if ~exist('cfg','var')
    cfg = [];
end
if ~exist('settings','var') | ~isfield(settings,'GRAPH') | ~isfield(settings.GRAPH,'measure')
    settings.GRAPH.folder = 'GraphResult';
    if exist('nummatrices','var') & ~isempty(nummatrices)
        settings.GRAPH.maxmatrix = nummatrices;
        if nummatrices > 1
            settings.GRAPH.avgmatrix = true;
            settings.GRAPH.numavg = nummatrices;
        else
            settings.GRAPH.avgmatrix = false;
            settings.GRAPH.numavg = [];
        end
    else
        settings.GRAPH.maxmatrix = [];
        settings.GRAPH.avgmatrix = false;
    end
    settings.GRAPH.Mappings = [];
    settings.GRAPH.MappingsMode = 'Average';
    settings.GRAPH.rankmatrix = false;
    settings.GRAPH.rankorder = 5;
    settings.GRAPH.measure = {};
    settings.GRAPH.MSTref = [];
    settings.GRAPH.MatrixRef = [];
    settings.GRAPH.randnumber = [];
    settings.GRAPH.randiter = [];
    if isfield(cfg,'Output_filepath') & ~isempty(cfg.Output_filepath)
        settings.GRAPH.deleteold = true;
    else
        settings.GRAPH.deleteold = false;
    end
end

Prompt = cell(0,2);
Formats = [];

Prompt(end+1,:) = {'Output-folder', 'folder'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'text';
Formats(end,1).size = 100;
Formats(end,1).span = [1 3];

if nonumbers == false
    Prompt(end+1,:) = {'Maximal number of matrices for graph analysis per subject (empty = all)','maxmatrix'};
    Formats(end+1,1).type = 'edit';
    Formats(end,1).format = 'integer';
    Formats(end,1).limits = [0 99999];
    Formats(end,1).size = 40;
    Formats(end,1).span = [1 3];
end

if ~exist('nummatrices','var') | isempty(nummatrices) | nummatrices > 1
    Prompt(end+1,:) = {'Average matrix per subject','avgmatrix'};
    Formats(end+1,1).type = 'check';
    
    Prompt(end+1,:) = {'Number of matrices for average (empty = all)','numavg'};
    Formats(end+1,1).type = 'edit';
    Formats(end,1).format = 'integer';
    Formats(end,1).limits = [0 99999];
    Formats(end,1).size = 40;
    Formats(end,1).span = [1 2];
end

Formats(end+1,1).type = 'none';
Formats(end,1).span = [1 3];

Prompt(end+1,:) = {'Reduce to Mappings','Mappings'};
Formats(end+1,1).type = 'check';
Formats(end,1).callback = {@load_mappings,'@ALL','@ALL',cfg};

Prompt(end+1,:) = {'Method','MappingsMode'};
Formats(end+1,1).type = 'list';
Formats(end,1).style = 'popupmenu';
Formats(end,1).format = 'input';
Formats(end,1).items = {'Average','Degree','Betweenness'};
Formats(end,1).span = [1 2];

Prompt(end+1,:) = {'Exclude channels','exclude'};
Formats(end+1,1).type = 'check';
Formats(end,1).callback = {@set_exclude,'@ALL','@ALL',header};
Formats(end,1).span = [1 3];

Prompt(end+1,:) = {'Rank matrix','rankmatrix'};
Formats(end+1,1).type = 'check';

Prompt(end+1,:) = {'Order for ranking','rankorder'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'float';
Formats(end,1).limits = [0 99];
Formats(end,1).size = 40;

Prompt{end+1,1} = {'',''};
Formats(end+1,1).type = 'text';

Formats(end+1,1).type = 'none';
Formats(end,1).span = [1 3];

Prompt(end+1,:) = {'Graph measures','measure'};
Formats(end+1,1).type = 'list';
Formats(end,1).style = 'listbox';
Formats(end,1).format = 'input';
Formats(end,1).items = {'ClusteringCoeff';'ClusteringCoeff_Norm';'ShortestPath'; ...
        'ShortestPath_Norm';'Degree';'DegreeCorrelation';'DegreeDiversity';'Density'; ...
        'Eccentricity';'Betweenness';'Efficiency';'EigenvectorCentrality'; ...
        'Transitivity';'Modularity';'Matrix-Ref'; ...
        'MST-Degree';'MST-DegreeCorrelation';'MST-DegreeDiversity';'MST-ShortestPath'; ...
        'MST-Eccentricity';'MST-Betweenness';'MST-Efficiency';'MST-EigenvectorCentrality'; ...
        'MST-LeaveNodes';'MST-Ref'};
Formats(end,1).limits = [0 2]; % multi-select
Formats(end,1).size = [170 355];
Formats(end,1).span = [15 2];
Formats(end,1).callback = {@lab_get_graphmeasures,'@ALL','@ALL'};

Prompt(end+1,:) = {'Number randomizations','randnumber'};
Formats(end+1,1).type = 'check';
Formats(end,1).callback = {@lab_get_number,'randnumber','randnumber'};

Prompt(end+1,:) = {'Iterations for randomizations','randiter'};
Formats(end+1,1).type = 'check';
Formats(end,1).callback = {@lab_get_number,'randiter','randiter'};

Prompt{end+1,1} = ' ';
Formats(end+1,1).type = 'text';

Prompt(end+1,:) = {'MST Ref','MSTref'};
Formats(end+1,1).type = 'check';
Formats(end,1).callback = {@lab_load_matrix,'MSTref','MSTref',[],true};

Prompt(end+1,:) = {'Matrix Ref','MatrixRef'};
Formats(end+1,1).type = 'check';
Formats(end,1).callback = {@lab_load_matrix,'MatrixRef','MatrixRef'};

Prompt{end+1,1} = ' ';
Formats(end+1,1).type = 'text';
Formats(end,1).span = [10 1];

Formats(end+1,1).type = 'none';
Formats(end,1).span = [1 3];

Prompt(end+1,:) = {'Delete results from previous run','deleteold'};
Formats(end+1,1).type = 'check';
Formats(end,1).span = [1 3];

[settings.GRAPH,Cancelled] = inputsdlg(Prompt,'Graph analysis',Formats,settings.GRAPH);
if isempty(settings.GRAPH) | Cancelled == 1
    skipprocessing = 1;
    settings.GRAPH = [];
end
pause(0.2);

end

function settings = load_mappings(settings,cfg)
    settings.Mappings = lab_load_mappings(settings.Mappings,cfg);
    if ~isempty(settings.Mappings)
        settings.exclude = [];
    end
end

function settings = set_exclude(settings,header)
    settings = lab_set_exclude(settings,header);
    if ~isempty(settings.exclude)
        settings.Mappings = [];
    end
end