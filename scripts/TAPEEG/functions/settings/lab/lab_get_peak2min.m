function [settings,skipprocessing] = lab_get_peak2min(settings,header,cfg)

skipprocessing = 0;

if ~exist('cfg','var')
    cfg = [];
end
if ~exist('header','var') | ~isfield(header,'numchannels')
    Channellist = cellstr(num2str((1:400)'));
    Locs = [];
else
    if isfield(header,'numdatachannels')
        includechannels = 1:header.numdatachannels;
    else
        includechannels = 1:header.numchannels;
    end
    if exist('cfg','var') & isfield(cfg,'exclude') & ~isempty(cfg.exclude) & ~isfield(header,'includechans')
        includechannels = setdiff(includechannels,cfg.exclude);
    end
    Channellist = cellstr(header.channels(includechannels,:));
    if isfield(header,'locs') & ~isempty(header.locs)
        Locs = header.locs;
        if exist('cfg','var') & isfield(cfg,'exclude') & ~isempty(cfg.exclude) & ~isfield(header,'includechans')
            Locs = lab_reduce_locs(Locs,cfg.exclude);
        end
    else
        Locs = [];
    end
end

if ~exist('settings','var') | ~isfield(settings,'BAchannels')
    settings.BAchannels = lab_get_bactivity(size(Channellist,1));
    settings.lowfreqpeak = 4;
    settings.highfreqpeak = 14;
    settings.MinPeak2Min = 1.3;
    settings.mode = 'threshold';
    settings.factor = 2;
    settings.threshold = 1.5;
end

Prompt = cell(0,2);
Formats = {};

Prompt(end+1,:) = {'Channels background-activity','BAchannels'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'vector';
Formats(end,1).limits = [-inf inf];
Formats(end,1).enable = 'inactive';
Formats(end,1).callback = {@lab_load_background,'BAchannels','BAchannels',Channellist,Locs};
Formats(end,1).span = [1 2];

Prompt(end+1,:) = {' ',''};
Formats(end+1,1).type = 'text';
Formats(end,1).span = [1 2];

Prompt(end+1,:) = {'Peak fequency: min', 'lowfreqpeak'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'float';
Formats(end,1).limits = [0 inf];
Formats(end,1).size = 40;

Prompt(end+1,:) = {'max', 'highfreqpeak'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'float';
Formats(end,1).limits = [0 inf];
Formats(end,1).size = 40;

Prompt(end+1,:) = {'Minimal Peak2Min for detection of peak frequency','MinPeak2Min'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'float';
Formats(end,1).limits = [0 inf];
Formats(end,1).size = 30;

Prompt(end+1,:) = {' ',''};
Formats(end+1,1).type = 'text';
Formats(end,1).span = [1 2];

Prompt(end+1,:) = {'Mode','mode'};
Formats(end+1,1).type = 'list';
Formats(end,1).style = 'popupmenu';
Formats(end,1).format = 'input';
Formats(end,1).items = {'threshold','weighted'};
Formats(end,1).callback = {@set_mode,'@ALL','@ALL'};
Formats(end,1).span = [1 2];

Prompt(end+1,:) = {'Threshold', 'threshold'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'float';
Formats(end,1).limits = [0 inf];
Formats(end,1).size = 40;

Prompt(end+1,:) = {'Reduction by', 'factor'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'float';
Formats(end,1).limits = [0 inf];
Formats(end,1).size = 40;
Formats(end,1).span = [1 2];

[settings,Cancelled] = inputsdlg(Prompt,'Detect by peak2min',Formats,settings);
if Cancelled == 1
    settings = [];
    skipprocessing = 1;
end

end

function settings = set_mode(settings)
   if strcmp(settings.mode,'weighted')
       settings.threshold = [];
       settings.factor = 3;
   elseif isempty(settings.threshold)
       settings.threshold = 1.5;
       settings.factor = 2;
   end
end