function [settings,skipprocessing] = lab_set_xls2matrix(settings,measures)

skipprocessing = 0;

if ~exist('measures','var')
    measures = cellstr(num2str((1:20)'));
end
if ~exist('settings','var') | ~isfield(settings,'folder')
    settings.folder = 'XLS2Matrix';
    settings.selection = 1;
    settings.invertvalues = 'off';
end

Prompt = cell(0,2);
Formats = {};

Prompt(end+1,:) = {'Output-folder', 'folder'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'text';
Formats(end,1).size = 100;
Formats(end,1).span = [1 2];

Prompt(end+1,:) = {'Measures','selection'};
Formats(end+1,1).type = 'list';
Formats(end,1).style = 'listbox';
Formats(end,1).items = measures;
Formats(end,1).limits = [0 2]; % multi-select
Formats(end,1).size = [140 160];
Formats(end,1).span = [1 2];

Prompt(end+1,:) = {'Invert values','invertvalues'};
Formats(end+1,1).type = 'list';
Formats(end,1).style = 'popupmenu';
Formats(end,1).format = 'input';
Formats(end,1).items = {'off','1-value','1/value'};
Formats(end,1).size = 140;

[settings,Cancelled] = inputsdlg(Prompt,'XLS to Matrix',Formats,settings);
if Cancelled == 1
    skipprocessing = 1;
end