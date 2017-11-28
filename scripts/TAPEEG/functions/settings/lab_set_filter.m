function [settings,skipprocessing] = lab_set_filter(settings,simple)

skipprocessing = 0;

if exist('simple','var') & isempty(simple)
    clearvars simple
end

if ~exist('settings','var') | isempty(settings) | ~isfield(settings,'filtermode')
    settings.filtermode = 'firls';
    settings.filtorder = 0;
    settings.wavsmoothing = [];
    if ~exist('simple','var')
        settings.notch = [];
    end
    settings.highpass = 0.5;
    settings.lowpass = 70;
    if ~exist('simple','var')
        settings.filterauxchannels = false;
    end
elseif strcmp(settings.filtermode,'freq') & isfield(settings,'dofreqshift') & settings.dofreqshift == 1
    settings.filtermode = 'freq-shift';
end
if isfield(settings,'off')
    settings = rmfield(settings,'off');
end

if ~exist('simple','var')
    disp ('   Ask for filter-settings')
end

Prompt = cell(0,2);
Formats = [];

Prompt(end+1,:) = {'Filter','filtermode'};
Formats(end+1,1).type = 'list';
Formats(end,1).style = 'popupmenu';
Formats(end,1).format = 'input';
Formats(end,1).items = {'firls','butter','freq','wavelet','cheby','freq-shift'};
Formats(end,1).callback = {@lab_get_filtermode,'@ALL','@ALL',1};

Prompt(end+1,:) = {'Order','filtorder'};
Formats(end+1,1).type = 'check';
Formats(end,1).callback = {@lab_get_filtermode,'@ALL','@ALL'};

Prompt(end+1,:) = {'WaveletSmoothing','wavsmoothing'};
Formats(end+1,1).type = 'check';
Formats(end,1).callback = {@lab_get_filtermode,'@ALL','@ALL'};

if ~exist('simple','var')
    Prompt(end+1,:) = {'Notch filter (Hz, empty = off)','notch'};
    Formats(end+1,1).type = 'edit';
    Formats(end,1).format = 'vector';
    Formats(end,1).limits = [0 inf];
    Formats(end,1).size = 60;
    Formats(end,1).span = [1 3];
end

Prompt(end+1,:) = {'Highpass filter (Hz, empty = off)','highpass'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'float';
Formats(end,1).limits = [0 inf];
Formats(end,1).size = 30;
Formats(end,1).span = [1 3];

Prompt(end+1,:) = {'Lowpass filter (Hz, empty = off)','lowpass'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'float';
Formats(end,1).limits = [0 inf];
Formats(end,1).size = 30;
Formats(end,1).span = [1 3];

if ~exist('simple','var')
    Prompt(end+1,:) = {'Filter aux channels','filterauxchannels'};
    Formats(end+1,1).type = 'check';
    Formats(end,1).size = [-1 -1];
    Formats(end,1).span = [1 3];
end

[settings,Cancelled] = inputsdlg(Prompt,'Filter',Formats,settings);
pause(0.2);
if isempty(settings) | Cancelled == 1
    settings = [];
    skipprocessing = 1;
else
    if strcmp(settings.filtermode,'freq-shift')
        settings.filtermode = 'freq';
        settings.dofreqshift = 1;
    else
        settings.dofreqshift = 0;
    end
end