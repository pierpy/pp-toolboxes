function [settings,skipprocessing] = lab_set_plot_xls2image(settings,Subjects,Measures,Vars,Vars2)

skipprocessing = 0;

settings.selectsubjects = 1:length(Subjects);
settings.selectvars = 1:length(Vars);
settings.selectvars2 = 1:length(Vars2);
settings.selectmeasures = 1:length(Measures);
settings.setsubjects = 1;
settings.setmeasures = 2;
settings.setvars = 3;
settings.setvars2 = 4;
settings.fontsize = 10;
settings.colorbar = false;
settings.flipcolormap = false;

if ~isempty(Vars) & ~isempty(Vars2)
    Items = {'Y-Axis','X-Axis','X-Blocks','Y-Blocks'};
elseif ~isempty(Vars)
    Items = {'Y-Axis','X-Axis','Blocks'};
else
    Items = {'Y-Axis','X-Axis'};
end

Prompt = cell(0,2);
Formats = {};

Prompt(end+1,:) = {'','setsubjects'};
Formats(end+1,1).type = 'list';
Formats(end,1).style = 'popupmenu';
Formats(end,1).items = Items;

Prompt(end+1,:) = {'','setmeasures'};
Formats(end+1,1).type = 'list';
Formats(end,1).style = 'popupmenu';
Formats(end,1).items = Items;

if ~isempty(Vars)
    Prompt(end+1,:) = {'','setvars'};
    Formats(end+1,1).type = 'list';
    Formats(end,1).style = 'popupmenu';
    Formats(end,1).items = Items;
end

if ~isempty(Vars2)
    Prompt(end+1,:) = {'','setvars2'};
    Formats(end+1,1).type = 'list';
    Formats(end,1).style = 'popupmenu';
    Formats(end,1).items = Items;
end

Prompt(end+1,:) = {'','selectsubjects'};
Formats(end+1,1).type = 'list';
Formats(end,1).style = 'listbox';
Formats(end,1).items = Subjects;
Formats(end,1).limits = [0 2]; % multi-select
Formats(end,1).size = [140 160];

Prompt(end+1,:) = {'','selectmeasures'};
Formats(end+1,1).type = 'list';
Formats(end,1).style = 'listbox';
Formats(end,1).items = Measures;
Formats(end,1).limits = [0 2]; % multi-select
Formats(end,1).size = [140 160];

if ~isempty(Vars)
    Prompt(end+1,:) = {'','selectvars'};
    Formats(end+1,1).type = 'list';
    Formats(end,1).style = 'listbox';
    Formats(end,1).items = Vars;
    Formats(end,1).limits = [0 2]; % multi-select
    Formats(end,1).size = [140 160];
else
    settings.selectvars = 1;
end

if ~isempty(Vars2)
    Prompt(end+1,:) = {'','selectvars2'};
    Formats(end+1,1).type = 'list';
    Formats(end,1).style = 'listbox';
    Formats(end,1).items = Vars2;
    Formats(end,1).limits = [0 2]; % multi-select
    Formats(end,1).size = [140 160];
else
    settings.selectvars2 = 1;
end

Prompt(end+1,:) = {'FontSize','fontsize'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'float';
Formats(end,1).limits = [0 inf];
Formats(end,1).size = 25;

Prompt(end+1,:) = {'Colorbar','colorbar'};
Formats(end+1,1).type = 'check';

Prompt(end+1,:) = {'Flip colormap','flipcolormap'};
Formats(end+1,1).type = 'check';

if ~isempty(Vars) & ~isempty(Vars2)
    [settings,Cancelled] = inputsdlg(Prompt,'Select data configuration',Formats,settings,4);
elseif ~isempty(Vars)
    [settings,Cancelled] = inputsdlg(Prompt,'Select data configuration',Formats,settings,3);
else
    [settings,Cancelled] = inputsdlg(Prompt,'Select data configuration',Formats,settings,2);
end
if Cancelled == 1
    skipprocessing = 1;
end