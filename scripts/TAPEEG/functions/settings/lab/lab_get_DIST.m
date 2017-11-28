function settings = lab_get_DIST(settings,simple)

if ~exist('simple','var')
    simple = false;
end

if ~exist('settings','var') | isempty(settings)
    settings.threshold = 5;
    settings.mode = '1';
    settings.input = 'off';
    settings.matrix = [];
end

Prompt = cell(0,2);
Formats = [];

if simple == false
    Prompt(end+1,:) = {'Matrix by name','input'};
    Formats(end+1,1).type = 'list';
    Formats(end,1).style = 'popupmenu';
    Formats(end,1).format = 'input';
    Formats(end,1).items = {'off','file','patient'};
    
    Prompt(end+1,:) = {'Default Matrix','matrix'};
    Formats(end+1,1).type = 'check';
    Formats(end,1).callback = {@lab_load_matrix,'matrix','matrix'};
else
    Prompt(end+1,:) = {'Distance Matrix','matrix'};
    Formats(end+1,1).type = 'check';
    Formats(end,1).callback = {@lab_load_matrix,'matrix','matrix'};
end

Prompt(end+1,:) = {'Threshold','threshold'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'float';
Formats(end,1).limits = [0 99999];
Formats(end,1).size = 30;

Prompt(end+1,:) = {'Value to set','mode'};
Formats(end+1,1).type = 'list';
Formats(end,1).style = 'popupmenu';
Formats(end,1).format = 'input';
Formats(end,1).items = {'1','max'};

[settings,Cancelled] = inputsdlg(Prompt,'Distance settings',Formats,settings);
if Cancelled == 1
    settings = [];
else
    pause(0.2);
end