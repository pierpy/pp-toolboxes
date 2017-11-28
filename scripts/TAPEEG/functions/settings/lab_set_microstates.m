function [cfg,skipprocessing] = lab_set_microstates(cfg,header)

skipprocessing = 0;

global THeader
if ~isempty(THeader)
    header = THeader;
elseif ~exist('header','var')
    header = [];
end

if ~exist('cfg','var') | isempty(cfg) | ~isfield(cfg,'MICROST') | ~isfield(cfg.MICROST,'eegsource')
    cfg.MICROST.folder = 'Microstates';
    cfg.MICROST.eegsource = 'input';
    cfg.MICROST.LAPL = [];
    cfg.MICROST.montage = [];
    cfg.MICROST.interpolate = true;
    cfg.MICROST.dogfp = true;
    if isfield(header,'samplingrate')
        cfg.MICROST.rejectsegm = round(header.samplingrate / 42);
    else
        cfg.MICROST.rejectsegm = 3;
    end
    cfg.MICROST.method = 'K-means';
    cfg.MICROST.minclusters = 2;
    cfg.MICROST.maxclusters = 20;
    cfg.MICROST.MaxIter = 300;
    cfg.MICROST.doplots = '2D';
    cfg.MICROST.doevents = [];
    cfg.MICROST.doevents_mode = 'fixed';
end

refitems = {'input','mean','median','laplacian','montage','channels'};
for i = 1:length(cfg.MICROST.eegsource)
    if ~max(strcmp({'input';'mean';'median';'laplacian';'montage'},cfg.MICROST.eegsource))
        refitems = cat(2,refitems,cellstr(cfg.MICROST.eegsource));
    end
end

Marker = lab_getall_markers(header,cfg);

Prompt = cell(0,3);
Formats = [];

Prompt(end+1,:) = {'Output-folder', 'folder',''};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'text';
Formats(end,1).size = 100;
Formats(end,1).span = [1 3];

Prompt(end+1,:) = {'Reference','eegsource',''};
Formats(end+1,1).type = 'list';
Formats(end,1).style = 'popupmenu';
Formats(end,1).format = 'input';
Formats(end,1).items = refitems;
Formats(end,1).callback = {@lab_get_eegsource,'@ALL','@ALL',cfg,header};

Prompt(end+1,:) = {'Montage','montage',''};
Formats(end+1,1).type = 'check';
Formats(end,1).callback = {@lab_load_montage,'montage','montage',cfg,header};

Prompt(end+1,:) = {'Laplacian','LAPL',''};
Formats(end+1,1).type = 'check';
Formats(end,1).callback = {@lab_get_laplacian,'LAPL','LAPL'};

Prompt(end+1,:) = {'Interpolate bad channels','interpolate',''};
Formats(end+1,1).type = 'check';
Formats(end,1).span = [1 3];

Formats(end+1,1).type = 'none';
Formats(end,1).span = [1 3];

Prompt{end+1,1} = 'Markers to exclude (strings / all)';
Formats(end+1,1).type = 'text';
Formats(end,1).span = [1 3];

Prompt(end+1,:) = {'', 'markerexclude',''};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'vector';
Formats(end,1).size = 300;
Formats(end,1).limits = [0 1];
if ~isempty(Marker)
    Formats(end,1).items = [{'all'} Marker(:)'];
end
Formats(end,1).callback = {@lab_get_marker,'markerexclude','markerexclude'};
Formats(end,1).span = [1 3];

Prompt(end+1,:) = {'Select by peaks in GFP','dogfp',''};
Formats(end+1,1).type = 'check';
Formats(end,1).span = [1 3];

Prompt(end+1,:) = {'Reject segments <','rejectsegm','TFs'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'integer';
Formats(end,1).limits = [1 9999];
Formats(end,1).size = 35;
Formats(end,1).span = [1 3];

Formats(end+1,1).type = 'none';
Formats(end,1).span = [1 3];

Prompt(end+1,:) = {'Template Clustering:','',''};
Formats(end+1,1).type = 'text';

Prompt(end+1,:) = {'Clustering:','',''};
Formats(end+1,1).type = 'text';

Prompt(end+1,:) = {'','method',''};
Formats(end+1,1).type = 'list';
Formats(end,1).style = 'popupmenu';
Formats(end,1).format = 'input';
Formats(end,1).items = {'K-means','T-AAHC'};
Formats(end,1).callback = {@set_method,'@ALL','@ALL'};

Prompt(end+1,:) = {'Template','template',''};
Formats(end+1,1).type = 'check';
Formats(end,1).callback = {@load_template,'@ALL','@ALL',header,cfg};

Prompt(end+1,:) = {'Collect all Files','doallfiles',''};
Formats(end+1,1).type = 'check';
Formats(end,1).format = 'integer';
Formats(end,1).callback = {@control_doallfiles,'@ALL','@ALL'};
Formats(end,1).span = [1 2];

Prompt(end+1,:) = {'Add template with zeros','AddZeros',''};
Formats(end+1,1).type = 'check';

Prompt(end+1,:) = {'Clusters','minclusters',''};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'integer';
Formats(end,1).limits = [2 9999];
Formats(end,1).size = 25;
Formats(end,1).callback = {@control_doevents,'@ALL','@ALL'};

Prompt(end+1,:) = {'-','maxclusters',''};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'integer';
Formats(end,1).limits = [2 9999];
Formats(end,1).size = 25;
Formats(end,1).callback = {@control_doevents,'@ALL','@ALL'};

Formats(end+1,1).type = 'none';

Prompt(end+1,:) = {'Maximal iterations','MaxIter',''};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'integer';
Formats(end,1).limits = [1 99999];
Formats(end,1).size = 35;
Formats(end,1).span = [1 2];

Formats(end+1,1).type = 'none';
Formats(end,1).span = [1 3];

Prompt(end+1,:) = {'Create Markers ','doevents','Clusters'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'integer';
Formats(end,1).limits = [2 9999];
Formats(end,1).size = 35;
Formats(end,1).callback = {@control_doevents,'@ALL','@ALL'};

Prompt(end+1,:) = {'','doevents_mode',''};
Formats(end+1,1).type = 'list';
Formats(end,1).style = 'popupmenu';
Formats(end,1).format = 'input';
Formats(end,1).items = {'disabled','fixed','auto first','auto'};
Formats(end,1).callback = {@control_doevents,'@ALL','@ALL'};
Formats(end,1).span = [1 2];

Formats(end+1,1).type = 'none';
Formats(end,1).span = [1 3];

Prompt(end+1,:) = {'Plot Clusters','doplots',''};
Formats(end+1,1).type = 'list';
Formats(end,1).style = 'popupmenu';
Formats(end,1).format = 'input';
Formats(end,1).items = {'disabled','2D','3D'};
Formats(end,1).span = [1 3];

[cfg.MICROST,Cancelled] = inputsdlg(Prompt,'Microstates',Formats,cfg.MICROST);
if isempty(cfg.MICROST) | Cancelled == 1
    cfg.MICROST = [];
    skipprocessing = 1;
else
    pause(0.2);
end

end

function settings = load_template(settings,header,cfg)
    settings.template = lab_load_microstates(settings.template,header);
    if isfield(header,'numdatachannels')
        Chans = 1:header.numdatachannels;
    elseif isfield(header,'numchannels')
        Chans = 1:header.numchannels;
    else
        Chans = [];
    end
    if isfield(cfg,'exclude') & ~isempty(cfg.exclude) & ~isempty(Chans)
        Chans = setdiff(Chans,cfg.exclude);
    end
    if ~isempty(Chans) & size(settings.template,1) ~= length(Chans)
        settings.template = [];
        warndlg('Number of channels in Microstates-File and Input-File not matching');
    end
    if ~isempty(settings.template)
        settings.maxclusters = [];
        settings.minclusters = [];
        settings.doevents = size(settings.template,2);
        settings.doevents_mode = 'fixed';
        settings.doallfiles = false;
    elseif isempty(settings.maxclusters)
        settings.minclusters = 2;
        settings.maxclusters = 20;
        settings.doevents = [];
        settings.doevents_mode = 'fixed';
        settings.doallfiles = false;
    end
end

function settings = control_doevents(settings)
    if ~isempty(settings.template)
        settings.minclusters = [];
        settings.maxclusters = [];
        settings.doevents = size(settings.template,2);
        settings.doevents_mode = 'fixed';
        settings.doallfiles = false;
        return
    end
    if settings.doevents > settings.maxclusters
        settings.doevents = settings.maxclusters;
    end
    if strcmp(settings.doevents_mode,'auto') | strcmp(settings.doevents_mode,'auto first')
        settings.doevents = [];
    end
end

function settings = control_doallfiles(settings)
    if settings.doallfiles == false
        settings.doallfiles = true;
        settings.doevents = [];
        settings.doevents_mode = 'disabled';
    else
        settings.doallfiles = false;
    end
    settings.template = [];
    if isempty(settings.minclusters)
        settings.minclusters = 2;
    end
    if isempty(settings.maxclusters)
        settings.maxclusters = 20;
    end
end

function settings = set_method(settings)
    if strcmp(settings.method,'T-AAHC')
        settings.MaxIter = [];
    elseif isempty(settings.MaxIter)
        settings.MaxIter = 300;
    end
end