function [cfg,skipprocessing] = lab_set_detect_eyeblinks(cfg,header)

skipprocessing = 0;

global THeader
if ~isempty(THeader)
    header = THeader;
elseif ~exist('header','var')
    header = [];
end

if ~exist('cfg','var')
    cfg = [];
end

if isfield(cfg,'EXTRA') & isfield(cfg.EXTRA,'numdatachans') & cfg.EXTRA.numdatachans == 257
    numchans = cfg.EXTRA.numdatachans;
elseif isfield(header,'numdatachannels')
    numchans = header.numdatachannels;
elseif isfield(header,'numchannels')
    numchans = header.numchannels;
else
    numchans = [];
end
if ~isfield(cfg,'EYEBLINKS') | ~isfield(cfg.EYEBLINKS,'eog') | isempty(cfg.EYEBLINKS.eog)
    if numchans == 257
        cfg.EYEBLINKS.eog = [32,257,0,0;37,257,0,0;46,257,0,0; ...
            25,257,0,0;18,257,0,0;10,257,0,0];
    elseif numchans == 214
        cfg.EYEBLINKS.eog = [31,214,0,0;36,214,0,0;45,214,0,0; ...
            25,214,0,0;18,214,0,0;10,214,0,0];
    elseif numchans == 204
        cfg.EYEBLINKS.eog = [31,204,0,0;36,204,0,0;45,204,0,0; ...
            25,204,0,0;18,204,0,0;10,204,0,0];
    else
        cfg.EYEBLINKS.eog = [];
    end
end
if ~isfield(cfg.EYEBLINKS,'SD') | isempty(cfg.EYEBLINKS.SD)
    cfg.EYEBLINKS.SD = 2;
end
if ~isfield(cfg.EYEBLINKS,'MinChans') | isempty(cfg.EYEBLINKS.MinChans)
    cfg.EYEBLINKS.MinChans = 66;
end
if ~isfield(cfg.EYEBLINKS,'duration') | isempty(cfg.EYEBLINKS.duration)
    cfg.EYEBLINKS.duration = 0.8;
end

Prompt = cell(0,2);
Formats = [];

Prompt(end+1,:) = {'Channel-list','eog'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'result';
Formats(end,1).callback = {@lab_get_eog,'eog','eog',header,cfg};
Formats(end,1).size = 70;

Formats(end+1,1).type = 'none';

Prompt(end+1,:) = {'Z-value','SD'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'float';
Formats(end,1).limits = [0 10];
Formats(end,1).size = 40;

Formats(end+1,1).type = 'none';

Prompt(end+1,:) = {'Minimum number of channels (%)','MinChans'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'float';
Formats(end,1).limits = [0 100];
Formats(end,1).size = 40;

Formats(end+1,1).type = 'none';

Prompt(end+1,:) = {'Minimal Duration (seconds)','duration'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'float';
Formats(end,1).limits = [0 10];
Formats(end,1).size = 40;

[cfg.EYEBLINKS,Cancelled] = inputsdlg(Prompt,'Detect Eyeblinks',Formats,cfg.EYEBLINKS);
if Cancelled == 1
    cfg.EYEBLINKS = [];
    return
end
pause(0.2);