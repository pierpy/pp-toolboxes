function settings = lab_get_SL(settings,dofreqs)

if ~exist('dofreqs','var')
    dofreqs = 0;
end
if exist('settings','var') & isfield(settings,'measure') & ~any(strcmp(settings.measure,'SL'))
    settings.SL = [];
    return
end

% SL settings
if ~exist('settings','var') | ~isfield(settings,'SL') | ~isfield(settings.SL,'pref')
    settings.SL.pref = 0.01;
    if dofreqs == 1
        settings.SL.freqs = [8 13];
    end
    settings.SL.storesingle = true;
    settings.SL.storehits = false;
    settings.SL.step = 0;
    settings.SL.stepunit = 'phases of max frequency';
end

Prompt = cell(0,2);
Formats = [];

Prompt(end+1,:) = {'Pref (SL)','pref'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'float';
Formats(end,1).limits = [0 inf];
Formats(end,1).size = 30;
Formats(end,1).span = [1 2];

Prompt(end+1,:) = {'Average matrix','storesingle'};
Formats(end+1,1).type = 'check';
Formats(end,1).size = [-1 -1];
Formats(end,1).span = [1 2];

Prompt(end+1,:) = {'Store hits','storehits'};
Formats(end+1,1).type = 'check';
Formats(end,1).size = [-1 -1];
Formats(end,1).span = [1 2];

if dofreqs == 1
    Prompt(end+1,:) = {'Freq range','freqs'};
    Formats(end+1,1).type = 'edit';
    Formats(end,1).format = 'vector';
    Formats(end,1).limits = [0 inf];
    Formats(end,1).size = 30;
    Formats(end,1).span = [1 2];
end

Prompt(end+1,:) = {'Stepsize (0 = no steps)','step'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'float';
Formats(end,1).limits = [0 9999];
Formats(end,1).size = 40;
Formats(end,1).callback = {@set_step,'storesingle','step','storesingle'};

Prompt(end+1,:) = {'','stepunit'};
Formats(end+1,1).type = 'list';
Formats(end,1).style = 'popupmenu';
Formats(end,1).format = 'input';
Formats(end,1).items = {'phases of max frequency';'TFs';'seconds'};

[settings.SL,Cancelled] = inputsdlg(Prompt,'SL settings',Formats,settings.SL);
if isempty(settings.SL) | Cancelled == 1
    settings.SL = [];
    return
end

end

function storesingle = set_step(step,storesingle)
   if step > 0
       storesingle = true;
   end
end

