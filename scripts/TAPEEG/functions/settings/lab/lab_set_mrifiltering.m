function settings = lab_set_mrifiltering(settings)

if ~exist('settings','var') | isempty(settings)
    settings.ts = 0.0625;
    settings.iter = 5;
    settings.cond = 3;
end

Prompt = cell(0,2);
Formats = [];

Prompt(end+1,:) = {'Sampling rate', 'ts'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'float';
Formats(end,1).limits = [0 10];
Formats(end,1).size = 60;

Prompt(end+1,:) = {'Number of iterations', 'iter'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'float';
Formats(end,1).limits = [0 99];
Formats(end,1).size = 30;

Prompt(end+1,:) = {'Diffusion parameter', 'cond'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'float';
Formats(end,1).limits = [0 99];
Formats(end,1).size = 30;

[settings,Cancelled] = inputsdlg(Prompt,'MRI filter',Formats,settings);
if Cancelled == 1
    settings = [];
    return
end

end