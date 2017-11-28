function settings = lab_set_plot_histogram(settings,data,measures)

if ~exist('measures','var')
    measures = [];
end
if ~exist('data','var')
    data = [];
end
if ~exist('settings','var') | ~isfield(settings,'combine')
    settings.combine = false;
end
if ~isfield(settings,'excludezeros')
    settings.excludezeros = false;
end
if ~isfield(settings,'minval')
    tmp = data(:);
    tmp = tmp(~isnan(tmp));
    settings.minval = round(min(tmp*100))/100;
    settings.maxval = round(max(tmp*100))/100;
    clearvars tmp
end
if ~isfield(settings,'nbins')
    settings.nbins = 10;
end
if settings.combine == true
    settings.numvertical = 1;
else
    if ~isfield(settings,'numvertical') | isempty(settings.numvertical)
        if size(data,1) > 3
            settings.numvertical= ceil(size(data,1)/4);
        else
            settings.numvertical = 1;
        end
    end
end

Prompt = cell(0,2);
Formats = [];

if ~isempty(measures)
    settings.selection = 1:size(measures,1);
    Prompt(end+1,:) = {'Measures','selection'};
    Formats(end+1,1).type = 'list';
    Formats(end,1).style = 'listbox';
    Formats(end,1).limits = [0 2]; % multi-select
    Formats(end,1).size = [140 160];
    Formats(end,1).items = measures;
    Formats(end,1).span = [1 2];
    Formats(end,1).callback = {@set_measures,'@ALL','@ALL',data};
end

Prompt(end+1,:) = {'Combine measures' 'combine'};
Formats(end+1,1).type = 'check';
Formats(end,1).format = 'integer';
Formats(end,1).enable = 'inactive';
Formats(end,1).span = [1 2];
Formats(end,1).callback = {@set_combine,'@ALL','@ALL',data};

Prompt(end+1,:) = {'Exclude zeros' 'excludezeros'};
Formats(end+1,1).type = 'check';
Formats(end,1).span = [1 2];

Prompt(end+1,:) = {'Minimal value', 'minval'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'float';
Formats(end,1).limits = [-inf inf];
Formats(end,1).size = 30;

Prompt(end+1,:) = {'Maximal value', 'maxval'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'float';
Formats(end,1).limits = [-inf inf];
Formats(end,1).size = 30;

if size(data,1) > 1
    Prompt(end+1,:) = {'Vertical number of plots', 'numvertical'};
    Formats(end+1,1).type = 'edit';
    Formats(end,1).format = 'float';
    Formats(end,1).limits = [1 20];
    Formats(end,1).size = 30;
    Formats(end,1).span = [1 2];
end

Prompt(end+1,:) = {'Number of bars', 'nbins'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'integer';
Formats(end,1).limits = [1 9999];
Formats(end,1).size = 30;
Formats(end,1).span = [1 2];

[settings,Cancelled] = inputsdlg(Prompt,'Plot histogram',Formats,settings);
if Cancelled == 1
    settings = [];
end

end

function settings = set_measures(settings,data)
   tmp = data(settings.selection,:);
   tmp = tmp(:);
   tmp = tmp(~isnan(tmp));
   settings.minval = round(min(tmp*100))/100;
   settings.maxval = round(max(tmp*100))/100;
   clearvars tmp
   if length(settings.selection) > 4
       settings.numvertical= ceil(length(settings.selection)/4);
   else
       settings.numvertical = 1;
   end
end

function settings = set_combine(settings,data)
   if settings.combine == false
       settings.combine = true;
       settings.numvertical = 1;
   else
       settings.combine = false;
       if size(data,1) > 4
           settings.numvertical= ceil(size(data,1)/4);
       else
           settings.numvertical = 1;
       end
   end
end