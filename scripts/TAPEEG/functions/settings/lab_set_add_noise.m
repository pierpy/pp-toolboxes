function [cfg,skipprocessing] = lab_set_add_noise(cfg,header)

skipprocessing = 0;

global THeader
if ~isempty(THeader)
    header = THeader;
elseif ~exist('header','var')
    header = [];
end

if ~exist('cfg','var') | ~isfield(cfg,'NOISE') | isempty(cfg.NOISE)
    cfg.NOISE.mode = 3;
    cfg.NOISE.dB = 10;
    cfg.NOISE.coeff = 0;
end

Prompt = cell(0,2);
Formats = [];

Prompt(end+1,:) = {'Mode','mode'};
Formats(end+1,1).type = 'list';
Formats(end,1).style = 'togglebutton';
Formats(end,1).items = {'AR-model';'Real data';'White noise';'Pink noise'};
Formats(end,1).callback = {@select_mode,'@ALL','@ALL'};

Prompt(end+1,:) = {'dB (noise level)','dB'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'float';
Formats(end,1).limits = [0 inf]; % 2-digits (positive #)
Formats(end,1).size = 25;

Prompt(end+1,:) = {'Coefficient (AR-model)','coeff'};
Formats(end+1,1).type = 'list';
Formats(end,1).style = 'popupmenu';
Formats(end,1).format = 'input';
Formats(end,1).items = [0 5 10 15 20 25 30];

Prompt(end+1,:) = {'Real data','data'};
Formats(end+1,1).type = 'check';
Formats(end,1).callback = {@select_data,'@ALL','@ALL'};

[cfg.NOISE,Cancelled] = inputsdlg(Prompt,'Add noise',Formats,cfg.NOISE);
if isempty(cfg) | Cancelled == 1
    cfg.NOISE = [];
    skipprocessing = 1;
end

    function settings = select_mode(settings)
        if settings.mode == 1
            if settings.coeff == 0
                settings.coeff = 10;
            end
        else
            settings.coeff = 0;
        end
        if settings.mode == 2
            settings = select_data(settings);
        else
            settings.data = [];
        end
    end
    
    function settings = select_data(settings)
        if settings.mode == 2
            if isfield(settings,'data')
                settings.data = lab_get_data(settings.data);
            else
                settings.data = lab_get_data;
            end
            if ~isempty(settings.data) & isfield(header,'numdatachannels') & size(settings.data,1) ~= header.numdatachannels
                warndlg('Number of channels not matching','Reject data')
                settings.data = [];
            end
            clearvars datatmp headertmp
        else
            settings.data = [];
        end
    end
    
end