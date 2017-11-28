function settings = lab_set_smoothing_auto(settings,flagxls)
    
if ~exist('settings','var') | ~isfield(settings,'maxdistance')
    settings.maxdistance = 10;
    settings.weightmaxdistance = 5;
    settings.steps = 0.25;
    settings.loops = 100;
    settings.omitted = 10;
    settings.dorelativ = false;
    if exist('flagxls','var') & flagxls == true
        settings.method = '3D';
    else
        settings.method = 'spherical';
    end
end
    
Prompt = cell(0,2);
Formats = {};

Prompt(end+1,:) = {'Maximal Distance (to calculate)','maxdistance'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'float';
Formats(end,1).limits = [0 inf];
Formats(end,1).size = 35;

Prompt(end+1,:) = {'Distance steps','steps'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'float';
Formats(end,1).limits = [0 inf];
Formats(end,1).size = 35;

Prompt(end+1,:) = {'Weight for maximal distance (percent)','weightmaxdistance'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'float';
Formats(end,1).limits = [0 inf];
Formats(end,1).size = 35;
Formats(end,1).span = [1 2];

Prompt(end+1,:) = {'Number of Loops','loops'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'integer';
Formats(end,1).limits = [0 inf];
Formats(end,1).size = 35;
Formats(end,1).span = [1 2];

Prompt(end+1,:) = {'Omitted channels per loop','omitted'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'integer';
Formats(end,1).limits = [0 inf];
Formats(end,1).size = 35;
Formats(end,1).span = [1 2];

Prompt(end+1,:) = {'Method for interpolation','method'};
Formats(end+1,1).type = 'list';
Formats(end,1).style = 'popupmenu';
Formats(end,1).format = 'input';
Formats(end,1).items = {'spherical','3D'};
Formats(end,1).span = [1 2];

if exist('flagxls','var') & flagxls == true
    Prompt(end+1,:) = {'Take relative data for calculation','dorelativ'};
    Formats(end+1,1).type = 'check';
end

[settings,Cancelled] = inputsdlg(Prompt,'Smoothing',Formats,settings);
if isempty(settings) | Cancelled == 1
    settings = [];
end

end