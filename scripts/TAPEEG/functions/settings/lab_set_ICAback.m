function [cfg,skipprocessing] = lab_set_ICAback(cfg,header,nopreproc,doshort)

skipprocessing = 0;

global THeader
if ~isempty(THeader)
    header = THeader;
end

if ~exist('doshort','var')
    doshort = 0;
end
if ~exist('nopreproc','var')
    nopreproc = 0;
end
if ~exist('header','var')
    header = [];
end
if isfield(header,'activations')
    activations = cellstr(num2str((1:header.activations)'));
elseif isfield(header,'numchannels')
    activations = cellstr(num2str((1:header.numchannels)'));
else
    activations = cellstr(num2str((1:400)'));
end

if ~exist('cfg','var') | ~isfield(cfg,'ICABACK') | ~isfield(cfg.ICABACK,'dobacktransform')
    if isfield(header,'badchans') | doshort == 1
        cfg.ICABACK.FILEBAD = true;
    else
        cfg.ICABACK.FILEBAD = false;
    end
    cfg.ICABACK.BAD = [];
    cfg.ICABACK.ACTIVATIONS = [];
    cfg.ICABACK.dobacktransform = true;
    cfg.ICABACK.foldername = 'ICAresult';
    cfg.ICABACK.IS = [];
    cfg.ICABACK.PREPROCESSING  = [];
end

Prompt = cell(0,2);
Formats = [];

if doshort == 0
    Prompt(end+1,:) = {'Bad activations:',''};
    Formats(end+1,1).type = 'text';
    
    Prompt(end+1,:) = {'File information (exclude.txt | .info)','FILEBAD'};
    Formats(end+1,1).type = 'check';
    Formats(end,1).format = 'integer';
    Formats(end,1).callback = {@set_file_bad,'@ALL','@ALL',header,activations};
    
    Prompt(end+1,:) = {'Detect bad','BAD'};
    Formats(end+1,1).type = 'check';
    Formats(end,1).callback = {@set_detect_bad,'@ALL','@ALL',header,cfg};
    
    Prompt(end+1,:) = {'Fixed bad','ACTIVATIONS'};
    Formats(end+1,1).type = 'edit';
    Formats(end,1).format = 'vector';
    Formats(end,1).limits = [-inf inf];
    Formats(end,1).items = activations;
    Formats(end,1).callback = {@set_select_activations,'@ALL','@ALL'};
    
    Formats(end+1,1).type = 'none';
    
    Prompt(end+1,:) = {'ICA backtransform',''};
    Formats(end+1,1).type = 'text';
    
    Prompt(end+1,:) = {'Enable','dobacktransform'};
    Formats(end+1,1).type = 'check';
end

Prompt(end+1,:) = {'Output Folder','foldername'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'text';
Formats(end,1).size = 100;

Prompt(end+1,:) = {'Dipol Fitting (activations)','IS'};
Formats(end+1,1).type = 'check';
Formats(end,1).callback = {@lab_set_dipolfit,'@ALL','@ALL'};

if nopreproc == 0
    Prompt(end+1,:) = {'Preprocessing','PREPROCESSING'};
    Formats(end+1,1).type = 'check';
    Formats(end,1).callback = {@set_preprocessing,'PREPROCESSING','PREPROCESSING',header,cfg};
end

[cfg.ICABACK,Cancelled] = inputsdlg(Prompt,'ICA-Backtransform',Formats,cfg.ICABACK);
if isempty(cfg.ICABACK) | Cancelled == 1
    cfg.ICABACK = [];
    skipprocessing = 1;
    return
end

end

function settings = set_detect_bad(settings,header,cfg)
    settings.BAD = lab_set_detect_bad(settings.BAD,header,cfg,1,1,2,1,0,1,1,0,0,0);
    if ~isempty(settings.BAD)
        settings.ACTIVATIONS = [];
        settings.FILEBAD = false;
    end
end

function settings = set_preprocessing(settings,header,cfg)
    if isempty(settings)
        settings = cfg;
    end
    settings = lab_set_preprocessing(settings,header,0,1);
end

function settings = set_file_bad(settings,header,activations)
    if isfield(header,'badchans') & ~isempty(header.badchans)
        selection = header.badchans;
    else
        selection = [];
    end
    selection = listdlg('PromptString','Bad Activations','SelectionMode','multiple','Name','List', ...
          'ListString',activations,'InitialValue',selection,'CancelString','None','ListSize',[200 260]);
    if isfield(header,'badchans') & length(selection) == length(header.badchans) & ...
            min(header.badchans == selection) == 1
        settings.FILEBAD = true;
        settings.ACTIVATIONS = [];
        settings.BAD = [];
    elseif ~isempty(selection)
        settings.FILEBAD = false;
        settings.ACTIVATIONS = selection;
        settings.BAD = [];
    else
        settings.FILEBAD = false;
        settings.ACTIVATIONS = [];
        settings.BAD = [];
    end
end

function settings = set_select_activations(settings)
    if ~isempty(settings.ACTIVATIONS)
        settings.BAD = [];
        settings.FILEBAD = false;
    end
end
