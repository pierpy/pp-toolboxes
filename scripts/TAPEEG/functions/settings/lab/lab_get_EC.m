function settings = lab_get_EC(settings)

if exist('settings','var') & isfield(settings,'measure') & ~any(strcmp(settings.measure,'EC'))
    settings.EC = [];
    return
end
    
if ~exist('settings','var') | ~isfield(settings,'EC') | ~isfield(settings.EC,'stepec')
    settings.EC.stepec = 0.5;
    settings.EC.EChilbert = true;
end

Prompt = cell(0,2);
Formats = [];

Prompt(end+1,:) = {'Step in seconds','stepec'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'float';
Formats(end,1).limits = [0 inf]; % 2-digits (positive #)
Formats(end,1).size = 30;

Prompt(end+1,:) = {'Use hilbert transform','EChilbert'};
Formats(end+1,1).type = 'check';
Formats(end,1).size = [-1 -1];

[settings.EC,Cancelled] = inputsdlg(Prompt,'EC settings',Formats,settings.EC);
if isempty(settings.EC) | Cancelled == 1
    settings.EC = [];
    return
end