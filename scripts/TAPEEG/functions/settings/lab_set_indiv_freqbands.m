function [cfg,skipprocessing] = lab_set_indiv_freqbands(cfg,header)

disp ('   Individual frequency bands settings')

skipprocessing = 0;

global THeader

if ~isempty(THeader)
    header = THeader;
elseif ~exist('header','var')
    header= [];
end
if ~exist('cfg','var')
    cfg = [];
end

if isfield(header,'channels')
    Channels = cellstr(header.channels);
else
    Channels = [];
end
if isfield(header,'locs')
    LOCS = header.locs;
else
    LOCS = [];
end

if ~exist('cfg','var') | ~isfield(cfg,'IFREQ') | ~isfield(cfg.IFREQ,'eegsource')
    cfg.IFREQ.lowfreq = 5;
    cfg.IFREQ.highfreq = 14;
    cfg.IFREQ.winsize = 4;
    if isfield(header,'numdatachannels')
        cfg.IFREQ.channels = lab_get_bactivity(header.numdatachannels);
    elseif isfield(header,'numchannels')
        cfg.IFREQ.channels = lab_get_bactivity(header.numchannels);
    else
        cfg.IFREQ.channels = [];
    end
end

Prompt = cell(0,2);
Formats = [];

Prompt(end+1,:) = {'Window length for fft (sec)','winsize'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'float';
Formats(end,1).limits = [0 inf];
Formats(end,1).size = 35;
Formats(end,1).span = [1 2];

Formats(end+1,1).type = 'none';
Formats(end,1).span = [1 2];

Prompt(end+1,:) = {'Band for peak frequency search',''};
Formats(end+1,1).type = 'text';
Formats(end,1).span = [1 2];

Prompt(end+1,:) = {'Lowest frequency','lowfreq'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'float';
Formats(end,1).limits = [0 inf];
Formats(end,1).size = 35;

Prompt(end+1,:) = {'Highest frequency','highfreq'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'float';
Formats(end,1).limits = [0 inf];
Formats(end,1).size = 35;

Formats(end+1,1).type = 'none';
Formats(end,1).span = [1 2];

Prompt(end+1,:) = {'Channels','channels'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'vector';
Formats(end,1).limits = [-inf inf];
Formats(end,1).enable = 'inactive';
Formats(end,1).size = 350;
Formats(end,1).span = [1 2];
Formats(end,1).callback = {@lab_load_background,'channels','channels',Channels,LOCS};

[cfg.IFREQ,Cancelled] = inputsdlg(Prompt,'Individual Frequency Bands',Formats,cfg.IFREQ,2);
pause(0.2);
if isempty(cfg.IFREQ) | Cancelled == 1
    cfg.IFREQ = [];
    skipprocessing = 1;
    return
end

end
