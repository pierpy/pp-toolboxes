function settings = lab_get_SLc(settings,dofreqs)

if ~exist('dofreqs','var')
    dofreqs = 0;
end
if exist('settings','var') & isfield(settings,'measure') & ~any(strcmp(settings.measure,'SLc'))
    settings.SLc = [];
    return
end

% SL settings
if ~exist('settings','var') | ~isfield(settings,'SLc') | ~isfield(settings.SLc,'pref')
    settings.SLc.pref = 0.05;
    settings.SLc.Mcorrect = 1;
    if dofreqs == 1
        settings.SLc.freqs = [8 13];
    end
    settings.SLc.storesingle = true;
    settings.SLc.storehits = false;
    settings.SLc.step = 0;
    settings.SLc.stepunit = 'phases of max frequency';
end

Prompt = cell(0,2);
Formats = [];

Prompt(end+1,:) = {'Pref (SL)','pref'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'float';
Formats(end,1).limits = [0 inf];
Formats(end,1).size = 30;
Formats(end,1).span = [1 2];

Prompt(end+1,:) = {'Mode','Mcorrect'};
Formats(end+1,1).type = 'list';
Formats(end,1).style = 'popupmenu';
Formats(end,1).items = {'uni-directional symmetrical','uni-directional asymmetrical','bi-directional'};

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

[settings.SLc,Cancelled] = inputsdlg(Prompt,'SLc settings',Formats,settings.SLc);
if isempty(settings.SLc) | Cancelled == 1
    settings.SLc = [];
    return
end

end

function storesingle = set_step(step,storesingle)
   if step > 0
       storesingle = true;
   end
end

