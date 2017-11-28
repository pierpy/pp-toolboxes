function [cfg,skipprocessing] = lab_set_create_modeldata(cfg)

skipprocessing = 0;

if ~exist('cfg','var') | ~isfield(cfg,'MODEL')
    cfg.MODEL.output_folder = '';
    cfg.MODEL.chans = 1;
    cfg.MODEL.fs = 500;
    cfg.MODEL.numtf = [10 11 12 13 14 15 16];
    cfg.MODEL.type = 'AR';
    cfg.MODEL.SET.coeff = 10;
    cfg.MODEL.loops = 1;
    cfg.MODEL.doconnect = false;
    cfg.MODEL.randmatrix = false;
    cfg.MODEL.matrixformat = 2;
    cfg.MODEL.lag = 5;
    cfg.MODEL.K = 1;
    cfg.MODEL.NOISE = [];
    cfg.MODEL.CONNECT = [];
    cfg.MODEL.GRAPH = [];
end

Prompt = cell(0,2);
Formats = [];

Prompt(end+1,:) = {'Output folder','output_folder'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'dir';
Formats(end,1).size = 250;
Formats(end,1).span = [1 2];

Prompt(end+1,:) = {'Number of channels','chans'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'integer';
Formats(end,1).limits = [0 999];
Formats(end,1).size = 35;
Formats(end,1).callback = {@correct_connections,'connections','connections','chans'};

Prompt(end+1,:) = {'Samplingrate','fs'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'float';
Formats(end,1).limits = [0 inf];
Formats(end,1).size = 45;

Prompt(end+1,:) = {'',''};
Formats(end+1,1).type = 'text';
Formats(end,1).span = [1 2];

Prompt(end+1,:) = {'Number of timeframes','numtf'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'vector';
Formats(end,1).limits = [0 inf];
Formats(end,1).size = 200;
Formats(end,1).span = [1 2];

Prompt(end+1,:) = {'',''};
Formats(end+1,1).type = 'text';
Formats(end,1).span = [1 2];

Prompt(end+1,:) = {'Model to use','type'};
Formats(end+1,1).type = 'list';
Formats(end,1).style = 'popupmenu';
Formats(end,1).format = 'input';
Formats(end,1).items = {'AR';'Freeman';'Sinuswave';'PhaseResetting';'Roessler'};
Formats(end,1).callback = {@lab_get_modelsettings,'@ALL','@ALL'};

Prompt(end+1,:) = {'Settings','SET'};
Formats(end+1,1).type = 'check';
Formats(end,1).callback = {@lab_get_modelsettings,'@ALL','@ALL'};

Prompt(end+1,:) = {'',''};
Formats(end+1,1).type = 'text';
Formats(end,1).span = [1 2];

Prompt(end+1,:) = {'Connections','connections'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'result';
Formats(end,1).callback = {@set_connections,'connections','connections','chans'};
Formats(end,1).size = 100;

Prompt(end+1,:) = {'Included','matrixformat'};
Formats(end+1,1).type = 'list';
Formats(end,1).style = 'popupmenu';
Formats(end,1).items = {'full','lower tirangle','upper triangle'};

Prompt(end+1,:) = {'Randomize connections','randmatrix'};
Formats(end+1,1).type = 'check';
Formats(end,1).size = [-1 -1];
Formats(end,1).span = [1 2];

Prompt(end+1,:) = {'Lag (timeframes / -1 = matrixvalues)','lag'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'integer';
Formats(end,1).limits = [-1 999];
Formats(end,1).size = 35;

Prompt(end+1,:) = {'K (coupling strength)','K'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'float';
Formats(end,1).limits = [0 1];
Formats(end,1).size = 35;

Prompt(end+1,:) = {'',''};
Formats(end+1,1).type = 'text';
Formats(end,1).span = [1 2];

Prompt(end+1,:) = {'Add noise','NOISE'};
Formats(end+1,1).type = 'check';
Formats(end,1).size = [-1 -1];
Formats(end,1).callback = {@lab_set_add_noise,'@ALL','@ALL'};
Formats(end,1).span = [1 2];

Prompt(end+1,:) = {'Calculate connectivity','CONNECT'};
Formats(end+1,1).type = 'check';
Formats(end,1).size = [-1 -1];
Formats(end,1).callback = {@lab_set_calculate_connectivity,'@ALL','@ALL'};

Prompt(end+1,:) = {'',''};
Formats(end+1,1).type = 'text';
Formats(end,1).span = [1 2];

Prompt(end+1,:) = {'Number of repetitions','loops'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'integer';
Formats(end,1).limits = [0 999];
Formats(end,1).size = 50;

[cfg.MODEL,Cancelled] = inputsdlg(Prompt,'Create model',Formats,cfg.MODEL);
if isempty(cfg.MODEL) | Cancelled == 1
    cfg.MODEL = [];
    skipprocessing = 1;
    return
else
    if min(cfg.MODEL.numtf) < 20
        cfg.MODEL.numtf = 2.^cfg.MODEL.numtf;
    end
end

end

function connections = correct_connections(connections,chans)
    if chans > 1
        if isempty(connections) | size(connections,2) ~= chans | size(connections,1) ~= chans
            connections = zeros(chans,chans);
        end
    else
        connections = [];
    end
end

function connections = set_connections(connections,nchans)
    if nchans > 1
        connections = lab_load_matrix(connections,nchans);
    else
        connections = [];
    end
end
 


