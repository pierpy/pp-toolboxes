function [cfg,skipprocessing] = lab_set_ICAstart(cfg,header,dobacktransform)

skipprocessing = 0;

global THeader
if ~isempty(THeader)
    header = THeader;
elseif ~exist('header','var')
    header = [];
end

if ~exist('dobacktransform','var')
    dobacktransform = true;
end

if ~exist('cfg','var') | ~isfield(cfg,'ICA') | ~isfield(cfg.ICA,'good')
    cfg.ICA.type = 'runica';
    cfg.ICA.extended = 0;
    cfg.ICA.pca = 0;
    cfg.ICA.good = true;
    cfg.ICA.epoch = 0;
    cfg.ICA.plottopo = false;
    cfg.ICA.BAD = [];
    cfg.ICA.BAD.length = 4;
    cfg.ICA.BAD.percentbad = 70;
    cfg.ICA.BAD.freqlim50 = 50;
    cfg.ICA.BAD.freqlim60 = 50;
    cfg.ICA.BAD.freqlimlow = 70;
    cfg.ICA.BAD.freqlimhigh = 50;
    cfg.ICA.BAD.spectslow = [0.5 2];
    cfg.ICA.BAD.spectshigh = [15 50];
    cfg.ICA.BAD.zvaluebroken = 0;
    cfg.ICA.BAD.zvaluevars = 0;
    cfg.ICA.BAD.zvaluehurst = 3;
    cfg.ICA.BAD.zvaluemedian = 3;
    cfg.ICA.BAD.zvaluekurtosis = 3;
    cfg.ICA.BAD.zvalueeye = 3;
    if exist('header','var') & isfield(header,'eog_ch') & header.eog_ch > 0
        cfg.ICA.BAD.eog = header.eog_ch;
    else
        if isfield(cfg,'EXTRA') & isfield(cfg.EXTRA,'numdatachans') & cfg.EXTRA.numdatachans == 257
            numchans = cfg.EXTRA.numdatachans;
        elseif isfield(header,'numdatachannels')
            numchans = header.numdatachannels;
        elseif isfield(header,'numchannels')
            numchans = header.numchannels;
        else
            numchans = [];
        end
        if numchans == 257
            cfg.ICA.BAD.eog = [37,241,46,244;18,238,10,234;31,252,32,253;31,226,25,225];
        elseif numchans == 214
            cfg.ICA.BAD.eog = [37,214,0,0;18,214,0,0];
        elseif numchans == 204
            cfg.ICA.BAD.eog = [37,204,0,0;18,204,0,0];
        else
            cfg.ICA.BAD.eog = [];
        end
    end
    if exist('header','var') & isfield(header,'ecg_ch') & ~isempty(header.ecg_ch) & header.ecg_ch ~= 0
        cfg.ICA.BAD.ecgdetect = 1;
        cfg.ICA.BAD.ecg_ch = header.ecg_ch;
    elseif exist('cfg','var') & isfield(cfg,'EXTRA') & isfield(cfg.EXTRA,'ecg_ch') & ...
            ~isempty(cfg.EXTRA.ecg_ch) & cfg.EXTRA.ecg_ch ~= 0
        cfg.ICA.BAD.ecgdetect = 1;
        cfg.ICA.BAD.ecg_ch = cfg.EXTRA.ecg_ch;
    else
        cfg.ICA.BAD.ecgdetect = 2;
        cfg.ICA.BAD.ecg_ch = [50 120];
    end
    cfg.ICA.BAD.PEAK2MIN = [];
    cfg.ICA.BAD.AVG = [];
    cfg.ICA.BAD.AVGstd = [];
    if isfield(cfg,'ICABACK')
        cfg.ICA.ICABACK = cfg.ICABACK;
    else
        cfg.ICA.ICABACK = [];
    end
end

Prompt = cell(0,2);
Formats = [];

Prompt(end+1,:) = {'Include only good channels','good'};
Formats(end+1,1).type = 'check';
Formats(end,1).size = [-1 -1];
Formats(end,1).span = [1 3];

Prompt(end+1,:) = {'Select epoch to process (0=all / samples / x-y)','epoch'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'vector';
Formats(end,1).limits = [-inf inf];
Formats(end,1).span = [1 3];

Formats(end+1,1).type = 'none';
Formats(end,1).span = [1 3];

Prompt(end+1,:) = {'Type','type'};
Formats(end+1,1).type = 'list';
Formats(end,1).style = 'popupmenu';
Formats(end,1).format = 'input'; % Answer will give value shown in items, disable to get integer
Formats(end,1).items = {'runica','fastICA','jader','binica'};

Prompt(end+1,:) = {'Extended (0 = off)','extended'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'integer';
Formats(end,1).limits = [-999 999];
Formats(end,1).size = 40;

Prompt(end+1,:) = {'PCA (0 = off, -1 = rank)','pca'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'integer';
Formats(end,1).limits = [0 999];
Formats(end,1).size = 40;

Formats(end+1,1).type = 'none';
Formats(end,1).span = [1 3];

Prompt(end+1,:) = {'Detect bad activations','BAD'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'result';
Formats(end,1).callback = {@lab_set_detect_bad,'BAD','BAD',header,cfg,1,1,2,1,0,1,1,0,0,0};
Formats(end,1).size = [150 280];
Formats(end,1).span = [1 3];

Prompt(end+1,:) = {'Plot ICA-Topo''s','plottopo'};
Formats(end+1,1).type = 'check';
Formats(end,1).span = [1 3];

Formats(end+1,1).type = 'none';
Formats(end,1).span = [1 3];

if dobacktransform == true
    Prompt(end+1,:) = {'Proceed automaticaly after detection of bad activations','ICABACK'};
    Formats(end+1,1).type = 'check';
    Formats(end,1).span = [1 3];
    Formats(end,1).callback = {@lab_set_ICAback,'@ALL','@ALL',header,1,1};
end

[cfg.ICA,Cancelled] = inputsdlg(Prompt,'ICA settings',Formats,cfg.ICA);
if isempty(cfg.ICA) | Cancelled == 1
    cfg.ICA = [];
    skipprocessing = 1;
    return
end

end
       
