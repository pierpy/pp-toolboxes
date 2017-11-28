function [settings,skipprocessing] = lab_set_eeg_Freeman(settings,mode)

skipprocessing = 0;

if ~exist('mode','var')
    mode = 0;
end

if ~exist('settings','var') | isempty(settings) | ~isfield(settings,'ae')
    settings = [];
    if mode == 1
        settings.chans = 1;
        settings.numtf = 4000;
        settings.fs = 1000;
    end
    settings.ae = 25;
    settings.be = 175;
    settings.ai = 75;
    settings.bi = 225;
    settings.vd = 15;
    settings.Qmax = 250;
    settings.kei = 5.0;
    settings.kie = 0.2;
    settings.sigma = 1;
    settings.tref = 5;
    settings.skip = 100;
end

Prompt = cell(0,2);
Formats = [];

if mode == 1
    Prompt(end+1,:) = {'Number channels','chans'};
    Formats(end+1,1).type = 'edit';
    Formats(end,1).format = 'integer';
    Formats(end,1).limits = [1 inf];
    Formats(end,1).size = 35;
    Formats(end,1).span = [1 2];
    
    Prompt(end+1,:) = {'Timeframes','numtf'};
    Formats(end+1,1).type = 'edit';
    Formats(end,1).format = 'integer';
    Formats(end,1).limits = [1 inf];
    Formats(end,1).size = 35;
    
    Prompt(end+1,:) = {'Samplingrate','fs'};
    Formats(end+1,1).type = 'edit';
    Formats(end,1).format = 'integer';
    Formats(end,1).limits = [1 inf];
    Formats(end,1).size = 35;
    
    Formats(end+1,1).type = 'none';
    Formats(end,1).span = [1 2];
end

Prompt(end+1,:) = {'Excitatory cell:',''};
Formats(end+1,1).type = 'text';
Formats(end,1).span = [1 2];

Prompt(end+1,:) = {'Rise time','ae'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'integer';
Formats(end,1).limits = [1 inf];
Formats(end,1).size = 35;

Prompt(end+1,:) = {'Decay time','be'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'integer';
Formats(end,1).limits = [1 inf];
Formats(end,1).size = 35;

Formats(end+1,1).type = 'none';
Formats(end,1).span = [1 2];

Prompt(end+1,:) = {'Excitatory cell:',''};
Formats(end+1,1).type = 'text';
Formats(end,1).span = [1 2];

Prompt(end+1,:) = {'Rise time','ai'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'integer';
Formats(end,1).limits = [1 inf];
Formats(end,1).size = 35;

Prompt(end+1,:) = {'Decay time','bi'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'integer';
Formats(end,1).limits = [1 inf];
Formats(end,1).size = 35;

Formats(end+1,1).type = 'none';
Formats(end,1).span = [1 2];

Prompt(end+1,:) = {'Threshold mV','vd'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'float';
Formats(end,1).limits = [0 inf];
Formats(end,1).size = 35;
Formats(end,1).span = [1 2];

Prompt(end+1,:) = {'Max firing rate','Qmax'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'float';
Formats(end,1).limits = [0 inf];
Formats(end,1).size = 35;
Formats(end,1).span = [1 2];

Formats(end+1,1).type = 'none';
Formats(end,1).span = [1 2];

Prompt(end+1,:) = {'Coupling:',''};
Formats(end+1,1).type = 'text';
Formats(end,1).span = [1 2];

Prompt(end+1,:) = {'E to I','kei'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'float';
Formats(end,1).limits = [0 inf];
Formats(end,1).size = 35;

Prompt(end+1,:) = {'I to E','kie'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'float';
Formats(end,1).limits = [0 inf];
Formats(end,1).size = 35;

Formats(end+1,1).type = 'none';
Formats(end,1).span = [1 2];

Prompt(end+1,:) = {'Sigma','sigma'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'float';
Formats(end,1).limits = [0 inf];
Formats(end,1).size = 35;
Formats(end,1).span = [1 2];

Prompt(end+1,:) = {'Refractory period','tref'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'float';
Formats(end,1).limits = [0 inf];
Formats(end,1).size = 35;
Formats(end,1).span = [1 2];

[settings,Cancelled] = inputsdlg(Prompt,'Neural mass model',Formats,settings,2);
if isempty(settings) | Cancelled == 1
    settings = [];
    skipprocessing = 1;
else
    pause(0.2);
end

end