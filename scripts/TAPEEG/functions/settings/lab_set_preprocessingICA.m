function [cfg,Prompt,Formats] = lab_set_preprocessingICA(cfg,header,Prompt,Formats)

dodiag = 1;

global THeader

if ~exist('Formats','var')
    Formats = [];
end
if ~exist('Prompt','var')
    Prompt = cell(0,2);
else
    dodiag = 0;
end
if isempty(THeader) & exist('header','var')
    THeader = header;
end
if ~exist('cfg','var') | ~isfield(cfg,'ICABACK')
    cfg.ICABACK.dobacktransform = true;
    cfg.ICABACK.BAD = [];
    cfg.ICABACK.PREPROCESSING  = [];
    cfg.ICABACK.ACTIVATIONS = [];
    cfg.ICABACK.foldername = '';
end

Prompt(end+1,:) = {'Preprocessing',[]};
Formats(end+1,1).type = 'text';
Formats(end,1).size = [-1 0];
Formats(end,1).span = [1 4]; % item is 1 field x 4 fields

Prompt(end+1,:) = {'ICA backtransform','ICABACK'};
Formats(end+1,1).type = 'check';
Formats(end,1).callback = {@lab_set_ICAback,'@ALL','@ALL'};

Prompt(end+1,:) = {'Stitch folder','STITCHALL'};
Formats(end+1,1).type = 'check';
Formats(end,1).callback = {@lab_set_stitchingall,'@ALL','@ALL'};

if dodiag == 1
    [cfg,Cancelled] = inputsdlg(Prompt,'Postprocessing',Formats,cfg);
    if isempty(cfg) | Cancelled == 1
        Prompt = 1;
        cfg = [];
        return
    else
        Prompt = 0;
    end
    
    THeader = [];
end

