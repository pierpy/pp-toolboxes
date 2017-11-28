function [settings,skipprocessing] = lab_set_interpolate_bad(settings,header)

disp ('   Ask for interpolate bad')
skipprocessing = 0;

if ~exist('settings','var') | ~isfield(settings,'method') | isempty(settings.method)
    settings.method = 'spherical';
end
if ~exist('header','var')
    header = [];
elseif isfield(header,'badchans')
    settings.badchannels = header.badchans;
else
    settings.badchannels = [];
end

Prompt = cell(0,2);
Formats = [];

if ~isempty(header)
    Prompt(end+1,:) = {'Bad channels','badchannels'};
    Formats(end+1,1).type = 'edit';
    Formats(end,1).format = 'vector';
    Formats(end,1).limits = [-inf inf];
    if isfield(header,'locs') & ~isempty(header.locs)
        Formats(end,1).enable = 'inactive';
        Formats(end,1).callback = {@get_badchannels,'badchannels','badchannels',header.locs};
    end
end

Prompt(end+1,:) = {'Method','method'};
Formats(end+1,1).type = 'list';
Formats(end,1).style = 'popupmenu';
Formats(end,1).format = 'input'; % Answer will give value shown in items, disable to get integer
Formats(end,1).items = {'spherical','3D'};

[settings,Cancelled] = inputsdlg(Prompt,'Interpolate bad',Formats,settings);
if isempty(settings) | Cancelled == 1
    settings = [];
    skipprocessing = 1;
    return
end

end

function badchannels = get_badchannels(badchannels,LOCS)
    settings.indexed = badchannels;
    settings.Color = [1 1 1];
    settings.ColorIdx = [1 0 0];
    settings.LOCS = LOCS;
    settings.Title = 'Bad channels';
    badchannels = lab_plot_locs(settings,1,0,0);
end