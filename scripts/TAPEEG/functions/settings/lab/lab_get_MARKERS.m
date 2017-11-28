function settings = lab_get_MARKERS(settings,header,cfg)

if ~exist('cfg','var')
    cfg = [];
end
if ~exist('header','var')
    header = [];
end
if ~exist('settings','var') | ~isfield(settings,'MARKER') | ~isfield(settings.MARKER,'markerexclude')
    if isfield(cfg,'BADELEC') & isfield(cfg.BADELEC,'markbad') & cfg.BADELEC.markbad == true
        settings.MARKER.markerexclude = cellstr('BAD');
    else
        settings.MARKER.markerexclude = {};
    end
end
if ~isfield(settings.MARKER,'markerinclude')
    settings.MARKER.markerinclude = {};
end
if ~isfield(settings.MARKER,'markerinclude_combine')
    settings.MARKER.markerinclude_combine = true;
end

Marker = lab_getall_markers(header,cfg);

if isfield(cfg,'MICROST') & isfield(cfg.MICROST,'doevents_mode') & ~strcmp(cfg.MICROST.doevents_mode,'disabled')
    docorrelation = true;
    if ~isfield(settings.MARKER,'mincorr')
        settings.MARKER.mincorr = 0.5;
    end
elseif isfield(header,'CORR') & ~isempty(header.CORR)
    docorrelation = true;
    if ~isfield(settings.MARKER,'mincorr')
        settings.MARKER.mincorr = 0.5;
    end
else
    docorrelation = false;
    settings.MARKER.mincorr = [];
end
    
Prompt = cell(0,2);
Formats = [];

Prompt{end+1,1} = 'Marker to exclude (strings / all)';
Formats(end+1,1).type = 'text';
Formats(end,1).span = [1 2];

Prompt(end+1,:) = {'', 'markerexclude'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'vector';
Formats(end,1).size = 300;
Formats(end,1).limits = [0 1];
if ~isempty(Marker)
    Formats(end,1).items = [{'all'} Marker(:)'];
end
Formats(end,1).callback = {@lab_get_marker,'markerexclude','markerexclude'};
Formats(end,1).span = [1 2];

if isfield(settings,'measure') & any(strcmp(settings.measure,'DTF'))
    Prompt{end+1,1} = 'Marker to include (strings / all / start)';
else
    Prompt{end+1,1} = 'Marker to include (strings / all)';
end
Formats(end+1,1).type = 'text';
Formats(end,1).span = [1 2];

Prompt(end+1,:) = {'', 'markerinclude'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'vector';
Formats(end,1).size = 300;
Formats(end,1).limits = [0 1];
if ~isempty(Marker)
    if isfield(settings,'measure') & any(strcmp(settings.measure,'DTF'))
        Formats(end,1).items = [{'all'} {'start'} Marker(:)'];
    else
        Formats(end,1).items = [{'all'} Marker(:)'];
    end
end
Formats(end,1).callback = {@lab_get_marker,'markerexclude','markerexclude'};

Prompt(end+1,:) = {'Combine','markerinclude_combine'};
Formats(end+1,1).type = 'check';

if docorrelation == true
    Prompt(end+1,:) = {'Minimal Cluster-correlation','mincorr'};
    Formats(end+1,1).type = 'edit';
    Formats(end,1).format = 'float';
    Formats(end,1).limits = [-1 1];
    Formats(end,1).size = 40;
end

[settings.MARKER,Cancelled] = inputsdlg(Prompt,'Markers',Formats,settings.MARKER);
if isempty(settings.MARKER) | Cancelled == 1
    settings.MARKER = [];
end

end