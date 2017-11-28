function [cfg,skipprocessing] = lab_set_locs(cfg,header)

skipprocessing = 0;

global THeader
if ~isempty(THeader)
    header = THeader;
elseif ~exist('header','var')
    header= [];
end

if ~exist('cfg','var') | ~isfield(cfg,'LOCS')
    cfg.LOCS = [];
end
if ~isfield(cfg.LOCS,'filelocs')
    if isfield(header,'locs') & ~isempty(header.locs)
        cfg.LOCS.filelocs = true;
    else
        cfg.LOCS.filelocs = false;
    end
end
if ~isfield(cfg.LOCS,'locs')
    cfg.LOCS.locs = [];
end
if ~isfield(cfg.LOCS,'matchlocs')
    cfg.LOCS.matchlocs = false;
end

% show dialog
if isfield(cfg,'MAIN') & isfield(cfg.MAIN,'auto') & cfg.MAIN.auto == 1
    return;
end

Prompt = cell(0,2);
Formats = [];

if isempty(header) | (isfield(header,'locs') & ~isfield(header.locs,'MATCHLOCS'))
    Prompt(end+1,:) = {'Use LOCS from eeg/meg-file','filelocs'};
    Formats(end+1,1).type = 'check';
    Formats(end,1).format = 'integer';
    Formats(end,1).callback = {@set_filelocs,'@ALL','@ALL'};
    Formats(end,1).span = [1 2];
    
    Formats(end+1,1).type = 'none';
    Formats(end+1,1).type = 'none';
end

Prompt(end+1,:) = {'Template LOCS-file','locs'};
Formats(end+1,1).type = 'check';
Formats(end,1).callback = {@load_locs,'@ALL','@ALL'};

Prompt(end+1,:) = {'Match LOCS','matchlocs'};
Formats(end+1,1).type = 'check';

[cfg.LOCS,Cancelled] = inputsdlg(Prompt,'Locs',Formats,cfg.LOCS,2);
if Cancelled == 1
    cfg.LOCS = [];
    skipprocessing = 1;
elseif isempty(cfg.LOCS.locs) & cfg.LOCS.filelocs == false
    cfg.LOCS = [];
end

end

function settings = load_locs(settings)
   settings.locs = lab_load_locs(settings.locs);
   if ~isempty(settings.locs)
       settings.filelocs = false;
   end
end

function settings = set_filelocs(settings)
   if settings.filelocs == false;
       settings.filelocs = true;
       settings.locs = [];
       settings.matchlocs = false;
   else
       settings.filelocs = false;
   end
end