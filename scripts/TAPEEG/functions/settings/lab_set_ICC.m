function [settings,skipprocessing] = lab_set_ICC(settings)

skipprocessing = 0;

if ~exist('settings','var') | isempty(settings) | ~isfield(settings,'mode')
    settings.type = 'agreement';
    settings.model = 'single-ratings';
    settings.mode = 'A-1';
    settings.orientation = 'measures';
    settings.alpha = 0.05;
    settings.ci = 'jackknife';
    settings.bootstrap = [];
    settings.plotdata = true;
end

Prompt = cell(0,2);
Formats = [];

Prompt(end+1,:) = {'ICC for','orientation'};
Formats(end+1,1).type = 'list';
Formats(end,1).style = 'popupmenu';
Formats(end,1).format = 'input';
Formats(end,1).items = {'measures','subjects'};

Formats(end+1,1).type = 'none';
Formats(end,1).span = [1 3];

Prompt(end+1,:) = {'Model','model'};
Formats(end+1,1).type = 'list';
Formats(end,1).style = 'popupmenu';
Formats(end,1).format = 'input';
Formats(end,1).items = {'single-ratings','mean-ratings'};
Formats(end,1).callback = {@set_mode,'@ALL','@ALL'};

Prompt(end+1,:) = {'Type','type'};
Formats(end+1,1).type = 'list';
Formats(end,1).style = 'popupmenu';
Formats(end,1).format = 'input';
Formats(end,1).items = {'agreement','consistency','one-way-random'};
Formats(end,1).callback = {@set_mode,'@ALL','@ALL'};

Prompt(end+1,:) = {'Mode','mode'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'text';
Formats(end,1).enable = 'inactive';
Formats(end,1).size = 60;

Formats(end+1,1).type = 'none';
Formats(end,1).span = [1 3];

Prompt(end+1,:) = {'Alpha','alpha'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'float';
Formats(end,1).limits = [0 1];
Formats(end,1).size = 30;
Formats(end,1).span = [1 3];

Formats(end+1,1).type = 'none';
Formats(end,1).span = [1 3];

Prompt(end+1,:) = {'Confidence interval','ci'};
Formats(end+1,1).type = 'list';
Formats(end,1).style = 'popupmenu';
Formats(end,1).format = 'input';
Formats(end,1).items = {'jackknife','bootstrap','off'};
Formats(end,1).callback = {@set_ci,'@ALL','@ALL'};

Prompt(end+1,:) = {'Loops','bootstrap'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'integer';
Formats(end,1).limits = [0 999999];
Formats(end,1).size = 50;
Formats(end,1).span = [1 2];

Formats(end+1,1).type = 'none';
Formats(end,1).span = [1 3];

Prompt(end+1,:) = {'Plot input data','plotdata'};
Formats(end+1,1).type = 'check';
Formats(end,1).span = [1 3];

[settings,Cancelled] = inputsdlg(Prompt,'ICC',Formats,settings);
if isempty(settings) | Cancelled == 1
    settings = [];
    skipprocessing = 1;
else
    pause(0.2);
end


    function settings = set_mode(settings)
        switch settings.type
            case 'consistency'
                settings.mode = 'C-';
            case 'agreement'
                settings.mode = 'A-';
            case 'one-way-random'
                settings.mode = '1-';
        end
        switch settings.model
            case 'single-ratings'
                settings.mode = [settings.mode '1'];
            case 'mean-ratings'
                settings.mode =  [settings.mode 'k'];
        end
    end

    function settings = set_ci(settings)
        switch settings.ci
            case 'jackknife'
                settings.bootstrap = [];
            case 'bootstrap'
                if isempty(settings.bootstrap)
                    settings.bootstrap = 10000;
                end
            case 'off'
                settings.bootstrap = [];
        end
    end
end