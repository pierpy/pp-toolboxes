function [settings,skipprocessing] = lab_set_CronbachAlpha(settings)

skipprocessing = 0;

if ~exist('settings','var') | isempty(settings) | ~isfield(settings,'ci')
    settings.orientation = 'measures';
    settings.ci = 'jackknife';
    settings.bootstrap = [];
    settings.plotdata = true;
end

Prompt = cell(0,2);
Formats = [];

Prompt(end+1,:) = {'Cronbach Alpha for','orientation'};
Formats(end+1,1).type = 'list';
Formats(end,1).style = 'popupmenu';
Formats(end,1).format = 'input';
Formats(end,1).items = {'measures','subjects'};
Formats(end,1).span = [1 2];

Formats(end+1,1).type = 'none';
Formats(end,1).span = [1 2];

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

Formats(end+1,1).type = 'none';
Formats(end,1).span = [1 2];

Prompt(end+1,:) = {'Plot input data','plotdata'};
Formats(end+1,1).type = 'check';
Formats(end,1).span = [1 2];

[settings,Cancelled] = inputsdlg(Prompt,'Cronbach Alpha',Formats,settings);
if isempty(settings) | Cancelled == 1
    settings = [];
    skipprocessing = 1;
else
    pause(0.2);
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