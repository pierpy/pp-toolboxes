function [cfg,Prompt,Formats] = lab_set_preprocessing(cfg,header,skipfilt,skipica,Prompt,Formats)

global THeader

dodiag = 1;
if ~exist('Prompt','var')
    Prompt = cell(0,2);
else
    dodiag = 0;
end
if ~exist('skipica','var')
    skipica = 0;
end
if ~exist('skipfilt','var')
    skipfilt = 0;
end
if ~exist('Formats','var')
    Formats = [];
end
if isempty(THeader) & exist('header','var')
    THeader = header;
end
if ~exist('cfg','var')
    cfg = [];
end

Prompt(end+1,:) = {'Preprocessing',[]};
Formats(end+1,1).type = 'text';
Formats(end,1).size = [-1 0];
Formats(end,1).span = [1 4]; % item is 1 field x 4 fields

if skipfilt == 0
    Prompt(end+1,:) = {'Scale','SCALE'};
    Formats(end+1,1).type = 'check';
    Formats(end,1).callback = {@lab_set_scale_data,'@ALL','@ALL'};
    
    Prompt(end+1,:) = {'Remove DBS','STIM'};
    Formats(end+1,1).type = 'check';
    Formats(end,1).callback = {@lab_set_remove_dbs,'@ALL','@ALL'};
    
    Prompt(end+1,:) = {'Filter','FILT'};
    Formats(end+1,1).type = 'check';
    Formats(end,1).callback = {@lab_set_filtering,'@ALL','@ALL'};
    
    Prompt(end+1,:) = {'Re-sample','RESAMPLE'};
    Formats(end+1,1).type = 'check';
    Formats(end,1).callback = {@lab_set_resample_data,'@ALL','@ALL'};
    
    if skipica == 0
        Prompt(end+1,:) = {'Stitch','STITCH'};
        Formats(end+1,1).type = 'check';
        Formats(end,1).callback = {@lab_set_stitching,'@ALL','@ALL'};
    end
    
    Prompt(end+1,:) = {'Detect bad','BADELEC'};
    Formats(end+1,1).type = 'check';
    Formats(end,1).callback = {@lab_set_find_good,'@ALL','@ALL'};
end

if skipica == 0
    Prompt(end+1,:) = {'ICA','ICA'};
    Formats(end+1,1).type = 'check';
    Formats(end,1).callback = {@lab_set_ICAstart,'@ALL','@ALL'};
end

Prompt(end+1,:) = {'Add noise','NOISE'};
Formats(end+1,1).type = 'check';
Formats(end,1).callback = {@lab_set_add_noise,'@ALL','@ALL'};

if skipica == 0
    Prompt(end+1,:) = {'Stitch folder','STITCHALL'};
    Formats(end+1,1).type = 'check';
    Formats(end,1).callback = {@lab_set_stitchingall,'@ALL','@ALL'};
end

Prompt(end+1,:) = {'Smoothing','SMOOTH'};
Formats(end+1,1).type = 'check';
Formats(end,1).callback = {@lab_set_smoothing,'@ALL','@ALL'};

Prompt(end+1,:) = {'Interpolate Bad','INTERPOLATE'};
Formats(end+1,1).type = 'check';
Formats(end,1).callback = {@lab_set_interpolate_bad,'INTERPOLATE','INTERPOLATE'};

Prompt(end+1,:) = {'Mark Eyeblinks','EYEBLINKS'};
Formats(end+1,1).type = 'check';
Formats(end,1).callback = {@lab_set_detect_eyeblinks,'@ALL','@ALL'};

Prompt(end+1,:) = {'Individual Frequency Bands','IFREQ'};
Formats(end+1,1).type = 'check';
Formats(end,1).callback = {@lab_set_indiv_freqbands,'@ALL','@ALL'};

if dodiag == 1
    [cfg,Cancelled] = inputsdlg(Prompt,'Preprocessing',Formats,cfg);
    if isempty(cfg) | Cancelled == 1
        cfg = [];
        Prompt = 1;
        return
    else
        Prompt = 0;
    end
    
    THeader = [];
end

