% commandline function for calculation of connectivity results
%
% [result,cfg] = lab_connect(data,header,cfg,nocfg)
%
% data & header     Output from lab_read_data
% cfg               config (if not set, config is created at first run)
% nocfg             if set & cfg is set from previous run, no dialog is shown
%
% written by F. Hatz 2013

function [Result,cfg] = lab_connect(data,header,cfg,nocfg)

Result = [];
if ~exist('header','var') | ~isfield(header,'samplingrate')
    header = lab_create_header(data);
end
if ~isfield(header,'numchannels')
    header.numchannels = size(data,1);
end
if ~isfield(header,'numdatachannels')
    header.numdatachannels = header.numchannels;
end
if ~isfield(header,'samplingrate')
    disp('Abort, samplingrate not defined')
end

if ~exist('cfg','var') | ~isfield(cfg,'CONNECT')
    cfg.CONNECT.filter='butter';
    cfg.CONNECT.spectralbands = lab_get_spectralbands;
    cfg.CONNECT.spectralbandsI = false;
    cfg.CONNECT.measure = {'PLI'};
    cfg.CONNECT.PHASE.step=0;
    cfg.CONNECT.PHASE.window= size(data,2) / header.samplingrate;
    cfg.CONNECT.PHASE.phaseestimate = 'hilbert';
    cfg.CONNECT.PHASE.usehann = true;
    cfg.CONNECT.PHASE.maxerror = 2;
    cfg.CONNECT.PHASE.delta = 2;
    cfg.CONNECT.PHASE.deltafreq = true;
    cfg.CONNECT.GRAPH = [];
    cfg.CONNECT.clustering = 0;
    cfg.CONNECT.correctdistance = 0;
    cfg.CONNECT.deleteold = true;
    [cfg,skipprocessing] = lab_set_calculate_connectivity(cfg,header,data,1);
    if skipprocessing == 1
        return
    end
elseif ~exist('nocfg','var')
    [cfg,skipprocessing] = lab_set_calculate_connectivity(cfg,header,data,1);
    if skipprocessing == 1
        return
    end
end
if isfield(cfg,'listold')
    cfg = rmfield(cfg,'listold');
end
[Result,cfg] = lab_calculate_connectivity(data,header,cfg);

return