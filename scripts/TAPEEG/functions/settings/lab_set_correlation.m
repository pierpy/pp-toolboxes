function [settings,skipprocessing] = lab_set_correlation(settings)

skipprocessing = 0;

if ~exist('settings','var') | isempty(settings) | ~isfield(settings,'type')
    settings.type = 'Pearson';
    settings.excludeNaN = false;
end

Prompt = cell(0,2);
Formats = {};
Prompt(end+1,:) = {'Select type','type'};
Formats(end+1,1).type = 'list';
Formats(end,1).style = 'popupmenu';
Formats(end,1).format = 'input';
Formats(end,1).items = {'Pearson','Kendall','Spearman'};

Formats(end+1,1).type = 'none';

Prompt(end+1,:) = {'Exclude NaN' 'excludeNaN'};
Formats(end+1,1).type = 'check';

[settings,Cancelled] = inputsdlg(Prompt,'Correlation',Formats,settings);
if Cancelled == 1
    skipprocessing = 1;
    return
end
