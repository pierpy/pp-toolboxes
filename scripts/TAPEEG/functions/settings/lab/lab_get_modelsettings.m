function settings = lab_get_modelsettings(settings)

Prompt = cell(0,2);
Formats = [];

if strcmp(settings.type,'AR')
    if ~isfield(settings,'SET') | ~isfield(settings.SET,'coeff') | isempty(settings.SET.coeff)
        settings.SET = [];
        settings.SET.coeff = 10;
    end
    Prompt(end+1,:) = {'Coefficient','coeff'};
    Formats(end+1,1).type = 'list';
    Formats(end,1).style = 'popupmenu';
    Formats(end,1).format = 'input';
    Formats(end,1).items = {5;10;15;20;25;30};    
elseif strcmp(settings.type,'Sinuswave')
    if ~isfield(settings,'SET') | ~isfield(settings.SET,'freq') | isempty(settings.SET.freq)
        settings.SET = [];
        settings.SET.freq = 8;
        settings.SET.randphase = true;
    end
    Prompt(end+1,:) = {'Frequency','freq'};
    Formats(end+1,1).type = 'edit';
    Formats(end,1).format = 'float';
    Formats(end,1).limits = [0 inf];
    Formats(end,1).size = 45;
    
    Prompt(end+1,:) = {'Randomize start phase','randphase'};
    Formats(end+1,1).type = 'check';
elseif strcmp(settings.type,'PhaseResetting')
    if ~isfield(settings,'SET') | ~isfield(settings.SET,'freq1') | isempty(settings.SET.freq1)
        settings.SET = [];
        settings.SET.sweepdur = 200;
        settings.SET.pfreq = 1;
        settings.SET.pamp = 11;
        settings.SET.ppos = 150;
        settings.SET.nfreq = 5;
        settings.SET.namp = 24;
        settings.SET.npos = 115;
        settings.SET.tjitter = 8;
        settings.SET.noiseamp = 10;
    end
    Prompt(end+1,:) = {'Duration (sweep)','sweepdur'};
    Formats(end+1,1).type = 'edit';
    Formats(end,1).format = 'integer';
    Formats(end,1).limits = [0 inf];
    Formats(end,1).size = 45;
    Formats(end,1).span = [1 3];
    
    Formats(end+1,1).type = 'none';
    Formats(end,1).span = [1 3];
    
    Prompt(end+1,:) = {'Positiv Peak: Frequency','pfreq'};
    Formats(end+1,1).type = 'edit';
    Formats(end,1).format = 'float';
    Formats(end,1).limits = [0 inf];
    Formats(end,1).size = 45;
    
    Prompt(end+1,:) = {'Amplitude','pamp'};
    Formats(end+1,1).type = 'edit';
    Formats(end,1).format = 'float';
    Formats(end,1).limits = [0 inf];
    Formats(end,1).size = 45;
    
    Prompt(end+1,:) = {'Position','ppos'};
    Formats(end+1,1).type = 'edit';
    Formats(end,1).format = 'integer';
    Formats(end,1).limits = [0 inf];
    Formats(end,1).size = 45;
    
    Prompt(end+1,:) = {'Negativ Peak: Frequency','nfreq'};
    Formats(end+1,1).type = 'edit';
    Formats(end,1).format = 'float';
    Formats(end,1).limits = [0 inf];
    Formats(end,1).size = 45;
    
    Prompt(end+1,:) = {'Amplitude','namp'};
    Formats(end+1,1).type = 'edit';
    Formats(end,1).format = 'float';
    Formats(end,1).limits = [0 inf];
    Formats(end,1).size = 45;
    
    Prompt(end+1,:) = {'Position','npos'};
    Formats(end+1,1).type = 'edit';
    Formats(end,1).format = 'integer';
    Formats(end,1).limits = [0 inf];
    Formats(end,1).size = 45;
    
    Formats(end+1,1).type = 'none';
    Formats(end,1).span = [1 3];
    
    Prompt(end+1,:) = {'Jitter','tjitter'};
    Formats(end+1,1).type = 'edit';
    Formats(end,1).format = 'integer';
    Formats(end,1).limits = [0 inf];
    Formats(end,1).size = 45;
    Formats(end,1).span = [1 3];
    
    Prompt(end+1,:) = {'Noise Amplitude','noiseamp'};
    Formats(end+1,1).type = 'edit';
    Formats(end,1).format = 'integer';
    Formats(end,1).limits = [0 inf];
    Formats(end,1).size = 45;
    Formats(end,1).span = [1 3];
elseif strcmp(settings.type,'Roessler')
    if ~isfield(settings,'SET') | ~isfield(settings.SET,'a') | isempty(settings.SET.a)
        settings.SET = [];
        settings.SET.a = 0.2;
        settings.SET.b = 0.4;
        settings.SET.c = 5.7;
        settings.SET.noise = 0;
    end
    Prompt(end+1,:) = {'a','a'};
    Formats(end+1,1).type = 'edit';
    Formats(end,1).format = 'float';
    Formats(end,1).limits = [0 inf];
    Formats(end,1).size = 45;
    
    Prompt(end+1,:) = {'b','b'};
    Formats(end+1,1).type = 'edit';
    Formats(end,1).format = 'float';
    Formats(end,1).limits = [0 inf];
    Formats(end,1).size = 45;
    
    Prompt(end+1,:) = {'c','c'};
    Formats(end+1,1).type = 'edit';
    Formats(end,1).format = 'float';
    Formats(end,1).limits = [0 inf];
    Formats(end,1).size = 45;
    
    Formats(end+1,1).type = 'none';
    
    Prompt(end+1,:) = {'Noise level','noise'};
    Formats(end+1,1).type = 'edit';
    Formats(end,1).format = 'float';
    Formats(end,1).limits = [0 inf];
    Formats(end,1).size = 45;
elseif strcmp(settings.type,'Freeman')
    settings.SET = lab_set_eeg_Freeman(settings.SET);
    return
end
if ~isempty(Prompt)
    [settings.SET,Cancelled] = inputsdlg(Prompt,[settings.type '-settings'],Formats,settings.SET);
    if isempty(settings.SET) | Cancelled == 1
        settings.SET = [];
    end
else
    settings.SET = [];
end
    
end