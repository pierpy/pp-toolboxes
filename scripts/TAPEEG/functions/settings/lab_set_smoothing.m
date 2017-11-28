function [cfg,skipprocessing] = lab_set_smoothing(cfg,header,flagxls)

skipprocessing = 0;

if ~exist('flagxls','var')
    flagxls = false;
end

global THeader
if ~isempty(THeader)
    header = THeader;
elseif ~exist('header','var')
    header = [];
end

if ~exist('cfg','var')
    cfg = [];
end

clearvars Prompt Formats
if ~isfield(cfg,'SMOOTH') | ~isfield(cfg.SMOOTH,'maxdistance')
    cfg.SMOOTH.folder = 'Smoothing';
    cfg.SMOOTH.Locs = [];
    cfg.SMOOTH.interpolate = true;
    cfg.SMOOTH.maxdistance = 4;
    cfg.SMOOTH.AUTO = [];
    cfg.SMOOTH.weightmaxdistance = 5;
end
if isfield(header,'locs') & ~isempty(header.locs)
    numchans = length(header.locs.x);
    skiplocs = true;
    cfg.SMOOTH.Locs = [];
elseif isfield(header,'channels')
    if isfield(header,'numdatachannels')
        numchans = header.numdatachannels;
    elseif isfield(header,'numhannels')
        numchans = header.numchannels;
    else
        numchans = size(header.channels,1);
    end
    skiplocs = false;
else
    numchans = [];
    skiplocs = false;
end

Prompt = cell(0,2);
Formats = {};

Prompt(end+1,:) = {'Output-folder', 'folder'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'text';
Formats(end,1).size = 100;
Formats(end,1).span = [1 4];

Formats(end+1,1).type = 'none';
Formats(end,1).span = [1 4];

if skiplocs == false
    Prompt(end+1,:) = {'LOCS','Locs'};
    Formats(end+1,1).type = 'check';
    Formats(end,1).callback = {@load_locs,'Locs','Locs'};
    Formats(end,1).span = [1 4];
    
    Formats(end+1,1).type = 'none';
    Formats(end,1).span = [1 4];
end

if flagxls == false
    Prompt(end+1,:) = {'Interpolate bad','interpolate'};
    Formats(end+1,1).type = 'check';
    Formats(end,1).span = [1 4];
end

Prompt(end+1,:) = {'Maximal Distance','maxdistance'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'float';
Formats(end,1).limits = [0 inf];
Formats(end,1).size = 35;
Formats(end,1).callback = {@set_maxdist,'@ALL','@ALL',flagxls};

Prompt(end+1,:) = {'Automatic','AUTO'};
Formats(end+1,1).type = 'check';
Formats(end,1).callback = {@set_auto,'@ALL','@ALL',flagxls};

Prompt(end+1,:) = {'Weight (percent)','weightmaxdistance'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'float';
Formats(end,1).limits = [0 inf];
Formats(end,1).size = 35;

[cfg.SMOOTH,Cancelled] = inputsdlg(Prompt,'Smoothing',Formats,cfg.SMOOTH);
if isempty(cfg.SMOOTH) | Cancelled == 1
    skipprocessing = 1;
    cfg.SMOOTH = [];
    return
else
    pause(0.2);
end
clearvars Prompt Formats

    function LOCS = load_locs(LOCS)
        LOCS = lab_load_locs(LOCS,cfg,numchans);
        if ~isempty(LOCS)
            numchans = length(LOCS.x);
        end
    end
    
end

function settings = set_auto(settings,flagxls)
   if ~isempty(settings.AUTO)
       settings.auto = false;
       if isempty(settings.maxdistance)
           settings.maxdistance = 4;
       end
   else
       settings.AUTO = lab_set_smoothing_auto(settings.AUTO,flagxls);
       if ~isempty(settings.AUTO)
           settings.maxdistance = [];
       end
   end
end

function settings = set_maxdist(settings,flagxls)
   if settings.maxdistance > 0
       settings.AUTO = [];
   else
       settings.AUTO = lab_set_smoothing_auto(settings.AUTO,flagxls);
       if ~isempty(settings.AUTO)
           settings.maxdistance = [];
       end
   end
end