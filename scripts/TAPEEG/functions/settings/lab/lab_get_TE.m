function settings = lab_get_TE(settings)

if exist('settings','var') & isfield(settings,'measure') & ~any(strcmp(settings.measure,'TE'))
    settings.TE = [];
    return
end

if ~exist('settings','var') | ~isfield(settings,'TE') | ~isfield(settings.TE,'el')
    settings.TE.el = 5;
    settings.TE.tau = 0.35;
    settings.TE.delay_min = 1;
    settings.TE.delay_step = 2;
    settings.TE.delay_max = 50;
end

Prompt = cell(0,2);
Formats = {};

Prompt(end+1,:) = {'Embedding length','el'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'integer';
Formats(end,1).limits = [0 inf];
Formats(end,1).size = 60;
Formats(end,1).span = [1 3];

Prompt(end+1,:) = {'Tau','tau'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'float';
Formats(end,1).limits = [0 inf];
Formats(end,1).size = 60;
Formats(end,1).span = [1 3];

Formats(end+1,1).type = 'none';
Formats(end,1).span = [1 3];

Prompt(end+1,:) = {'Delay',''};
Formats(end+1,1).type = 'text';
Formats(end,1).span = [1 3];

Prompt(end+1,:) = {'Min','delay_min'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'integer';
Formats(end,1).limits = [0 inf];
Formats(end,1).size = 40;

Prompt(end+1,:) = {'Step','delay_step'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'integer';
Formats(end,1).limits = [0 inf];
Formats(end,1).size = 40;

Prompt(end+1,:) = {'Max','delay_max'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'integer';
Formats(end,1).limits = [0 inf];
Formats(end,1).size = 40;

[settings.TE,Cancelled] = inputsdlg(Prompt,'TE',Formats,settings.TE);
if isempty(settings.TE) | Cancelled == 1
    settings.TE = [];
    return
end

end