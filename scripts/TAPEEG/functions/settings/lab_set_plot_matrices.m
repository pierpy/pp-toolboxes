function settings = lab_set_plot_matrices(settings,header,Mappings)

global Main_Path

if ~exist('header','var')
    header = [];
end
if ~exist('settings','var') | ~isfield(settings,'threshold')
    settings.mode = 'signal space';
    settings.LOCS = [];
    settings.ThresholdE = 96.88;
    settings.SizeE = 1;
    settings.ColorE = [0 0 1];
    settings.Size = 1;
    settings.Color = [1 0 0];
    settings.nodemethod = 'No nodes';
    settings.nodenormalize = false;
    settings.surfacemethod = 'No surface';
    settings.surfacenormalize = false;
    settings.MinValue = 0;
    settings.MaxValue = 0;
else
    settings.Size = settings.Size / 0.016;
    settings.Color = settings.Color(end,:);
end
if isempty(settings.LOCS)
    if isfield(header,'locs') & ~isempty(header.locs)
        settings.LOCS = header.locs;
    elseif exist(fullfile(Main_Path,'electrodes.els'),'file')
        settings.LOCS = lab_read_data(fullfile(Main_Path,'electrodes.els'));
    elseif exist(fullfile(Main_Path,'electrodes.sfp'),'file')
        settings.LOCS = lab_read_data(fullfile(Main_Path,'electrodes.sfp'));
    end
end

Prompt = cell(0,2);
Formats = [];

Prompt(end+1,:) = {'','mode'};
Formats(end+1,1).type = 'list';
Formats(end,1).style = 'popupmenu';
Formats(end,1).format = 'input';
Formats(end,1).items = {'signal space','source space'};
Formats(end,1).callback = {@set_mode,'@ALL','@ALL'};

Prompt(end+1,:) = {'Locs','LOCS'};
Formats(end+1,1).type = 'check';
Formats(end,1).callback = {@load_electrodes,'@ALL','@ALL'};

if exist('Mappings','var')
    settings.Mappings = Mappings;
    Prompt(end+1,:) = {'Mappings','Mappings'};
    Formats(end+1,1).type = 'check';
    Formats(end,1).size = [-1 -1];
    Formats(end,1).callback = {@lab_load_mappings,'Mappings','Mappings'};
end

Formats(end+1,1).type = 'none';
Formats(end,1).span = [1 3];

Prompt(end+1,:) = {'Threshold (%)','ThresholdE'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'float';
Formats(end,1).limits = [0 inf];
Formats(end,1).size = 50;
Formats(end,1).span = [1 3];

Prompt(end+1,:) = {'Edge size (default = 1)', 'SizeE'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'integer';
Formats(end,1).limits = [0 99];
Formats(end,1).size = 40;
Formats(end,1).span = [1 2];

Prompt(end+1,:) = {'Edge color','ColorE'};
Formats(end+1,1).type = 'color';
Formats(end,1).size = 20;

Formats(end+1,1).type = 'none';
Formats(end,1).span = [1 3];

Prompt(end+1,:) = {'Nodes','nodemethod'};
Formats(end+1,1).type = 'list';
Formats(end,1).style = 'popupmenu';
Formats(end,1).format = 'input';
Formats(end,1).items = {'No nodes','Degree','Betweenness centrality','Eigenvector centrality', ...
    'Clustering coefficient','Small dots','Large dots','Diagonal data'};
Formats(end,1).span = [1 2];

Prompt(end+1,:) = {'Normalize','nodenormalize'};
Formats(end+1,1).type = 'check';

Prompt(end+1,:) = {'Node size (default = 1)', 'Size'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'integer';
Formats(end,1).limits = [0 99];
Formats(end,1).size = 40;
Formats(end,1).span = [1 2];

Prompt(end+1,:) = {'Node color','Color'};
Formats(end+1,1).type = 'color';
Formats(end,1).size = 20;

Formats(end+1,1).type = 'none';
Formats(end,1).span = [1 3];

Prompt(end+1,:) = {'Min Value','MinValue'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'float';
Formats(end,1).limits = [0 inf];
Formats(end,1).size = 50;

Prompt(end+1,:) = {'Max Value','MaxValue'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'float';
Formats(end,1).limits = [0 inf];
Formats(end,1).size = 50;

Prompt(end+1,:) = {'(0 = automatic)',''};
Formats(end+1,1).type = 'text';

[settings,Cancelled] = inputsdlg(Prompt,'Plot matrix',Formats,settings);
if Cancelled == 1
    settings = [];
else
    settings.Size = settings.Size * 0.016;
    settings.add = 0;
end

end

function settings = load_electrodes(settings)
   if strcmp(settings.mode,'signal space')
       settings.LOCS = lab_load_locs(settings.LOCS);
   else
       settings.LOCS = [];
   end
end

function settings = set_mode(settings)
    if strcmp(settings.mode,'signal space')
        settings.ColorE = [0 0 1];
        settings.Color = [1 0 0];
    else
        settings.ColorE = [1 0 0];
        settings.Color = [0 0 0];
    end
end