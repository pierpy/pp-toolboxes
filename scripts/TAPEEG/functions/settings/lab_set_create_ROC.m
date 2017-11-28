function [settings,skipprocessing] = lab_set_create_ROC(settings,outcomes)

skipprocessing = 0;
if ~exist('outcomes','var') | isempty(outcomes)
    outcomes = logical([0 1]);
end

if ~exist('settings','var') | ~isfield(settings,'method')
    settings.method = '<=';
    settings.nbootstrap = 0;
    settings.jackknife = false;
end
if ~exist('settings','var') | ~isfield(settings,'target')
    if isnumeric(outcomes)
        List = unique(outcomes);
        List = num2cell(num2str(List(:)));
        settings.target = num2str(max(outcomes));
        doconvert = 1;
    elseif iscell(outcomes)
        List = unique(outcomes);
        List = List(:);
        settings.target = List{1};
        doconvert = 0;
    else
        settings.target = '1';
        List = {'1';'0'};
        doconvert = 2;
    end
elseif islogical(settings.target)
    settings.target = num2str(settings.target);
    List = {'1';'0'};
    doconvert = 2;
elseif isnumeric(settings.target)
    if settings.target > 0
        List = 1:settings.target;
        List = num2cell(num2str(List(:)));
        doconvert = 1;
    else
        List = {'1';'0'};
        doconvert = 2;
    end
    settings.target = num2str(settings.target);
elseif ischar(settings.target)
    doconvert = 0;
    List = cellstr(settings.target);
else
    settings = [];
    skipprocessing = 1;
    return
end

Prompt = cell(0,2);
Formats = [];

Prompt(end+1,:) = {'Rule for Positiv','method'};
Formats(end+1,1).type = 'list';
Formats(end,1).style = 'popupmenu';
Formats(end,1).format = 'input';
Formats(end,1).items = {'<=';'<';'>=';'>'};

Prompt(end+1,:) = {'Positive Outcome','target'};
Formats(end+1,1).type = 'list';
Formats(end,1).style = 'popupmenu';
Formats(end,1).format = 'input';
Formats(end,1).items = List(:);

Prompt(end+1,:) = {'Bootstrap (0 = off)','nbootstrap'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'integer';
Formats(end,1).limits = [0 999999999];
Formats(end,1).size = 80;

Prompt(end+1,:) = {'Jackknife' 'jackknife'};
Formats(end+1,1).type = 'check';

[settings,Cancelled] = inputsdlg(Prompt,'ROC-Curve',Formats,settings);
if Cancelled == 1
    settings = [];
    skipprocessing = 1;
elseif doconvert == 2
    settings.target = logical(str2num(settings.target));
elseif doconvert == 1
    settings.target = str2num(settings.target);
end

end