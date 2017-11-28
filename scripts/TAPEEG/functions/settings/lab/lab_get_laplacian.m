function [settings,skipprocessing] = lab_get_laplacian(settings,flag)

skipprocessing = 0;

clearvars Prompt Formats
if ~exist('settings','var') | ~isfield(settings,'lap_maxdistance')
    if exist('flag','var') & flag == 1
        settings.lap_maxdistance = 2.5;
        settings.lap_weightmaxdistance = 30;
        settings.lap_excluderef = true;
    else
        settings.lap_maxdistance = 4;
        settings.lap_weightmaxdistance = 5;
        settings.lap_excluderef = false;
    end
end
Prompt = cell(0,2);
Formats = {};
Prompt(end+1,:) = {'Maximal distance (1=minimal electrode distance)', 'lap_maxdistance'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'float';
Formats(end,1).limits = [0 inf];
Formats(end,1).size = 40;

Prompt(end+1,:) = {'Weight of maximal distance (percent)', 'lap_weightmaxdistance'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'float';
Formats(end,1).limits = [0 inf];
Formats(end,1).size = 40;

Prompt(end+1,:) = {'Exclude reference channel' 'lap_excluderef'};
Formats(end+1,1).type = 'check';
Formats(end,1).size = [-1 -1];

[settings,Cancelled] = inputsdlg(Prompt,'Laplacian reference',Formats,settings);
if Cancelled == 1
    settings = [];
    skipprocessing = 1;
end