function [cfg,skipprocessing] = lab_set_calculate_connectivity(cfg,header,data,nodel)

skipprocessing = 0;
disp ('   Ask for Connectivity-settings')

if ~exist('nodel','var')
    nodel = 0;
end

global THeader TData
if ~isempty(TData)
    data = TData;
elseif ~exist('data','var')
    data = [];
end
if ~isempty(THeader)
    header = THeader;
elseif ~exist('header','var')
    header = [];
end

if ~exist('cfg','var') | ~isfield(cfg,'CONNECT') | ~isfield(cfg.CONNECT,'filter')
    cfg.CONNECT.folder = 'Connectivity';
    cfg.CONNECT.interpolate = true;
    cfg.CONNECT.filter = 'butter';
    cfg.CONNECT.spectralbands = lab_get_spectralbands;
    cfg.CONNECT.spectralbands = cfg.CONNECT.spectralbands(2:end,:);
    cfg.CONNECT.spectralbandsI = false;
    cfg.CONNECT.measure = {'PLI'};
    cfg.CONNECT.doaverage = false;
    cfg.CONNECT.randphase = false;
    cfg.CONNECT.numrandphase = [];
    cfg.CONNECT.writematrix = true;
    cfg.CONNECT.PLOT = [];
    cfg.CONNECT.deleteold = true;
    cfg.CONNECT.correctdistance = 0;
    cfg.CONNECT.clustering = 0;
    cfg.CONNECT.PHASE.step = 4096;
    cfg.CONNECT.PHASE.stepunit = 'TFs';
    cfg.CONNECT.PHASE.phaseestimate = 'hilbert';
    cfg.CONNECT.PHASE.window = 4;
    cfg.CONNECT.PHASE.usehann = true;
    cfg.CONNECT.PHASE.debiased = true;
    cfg.CONNECT.PHASE.maxerror = 2;
    cfg.CONNECT.PHASE.flagold = false;
    cfg.CONNECT.PHASE.delta = 2;
    cfg.CONNECT.PHASE.deltafreq = true;
    cfg.CONNECT.PHASE.dosymm = false;
    cfg.CONNECT.PHASE.donorm = true;
    cfg.CONNECT.PHASE.downsample = 2;
end
if ~isfield(cfg.CONNECT,'eegsource') | isempty(cfg.CONNECT.eegsource)
    if isfield(header,'ref_chan') & ~isempty(header.ref_chan)
        if isnumeric(header.ref_chan) & ~isempty(header.ref_chan)
            cfg.CONNECT.eegsource = num2str(header.ref_chan);
        else
            cfg.CONNECT.eegsource = header.ref_chan;
        end
    else
        cfg.CONNECT.eegsource = 'input';
    end
end
if isfield(cfg.CONNECT,'spectralbands') & isnumeric(cfg.CONNECT.spectralbands) & size(cfg.CONNECT.spectralbands,2) == 2
    cfg.CONNECT.spectralbands = cat(2,cell(size(cfg.CONNECT.spectralbands,1),1),num2cell(cfg.CONNECT.spectralbands));
end

Prompt = cell(0,2);
Formats = [];

Prompt(end+1,:) = {'Output-folder', 'folder'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'text';
Formats(end,1).size = 100;
Formats(end,1).span = [1 4];

Prompt(end+1,:) = {'Reference','eegsource'};
Formats(end+1,1).type = 'list';
Formats(end,1).style = 'popupmenu';
Formats(end,1).format = 'input';
if ~max(strcmp({'channels','mean','median','laplacian','input','montage','input'},cfg.CONNECT.eegsource))
    Formats(end,1).items = {cfg.CONNECT.eegsource,'channels','mean','median','laplacian','montage','input'};
else
    Formats(end,1).items = {'channels','mean','median','laplacian','montage','input'};
end
Formats(end,1).callback = {@lab_get_eegsource,'@ALL','@ALL',cfg,header};

Prompt(end+1,:) = {'Montage','montage'};
Formats(end+1,1).type = 'check';
Formats(end,1).callback = {@lab_load_montage,'montage','montage',cfg,header};

Prompt(end+1,:) = {'Laplacian','LAPL'};
Formats(end+1,1).type = 'check';
Formats(end,1).callback = {@lab_get_laplacian,'LAPL','LAPL'};

Prompt(end+1,:) = {'Interpolate bad channels','interpolate'};
Formats(end+1,1).type = 'check';

Formats(end+1,1).type = 'none';
Formats(end,1).span = [1 4];

Prompt(end+1,:) = {'Spectral Bands','spectralbands'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'result';
Formats(end,1).callback = {@lab_get_spectralbands,'spectralbands','spectralbands','spectralbandsI',true};
Formats(end,1).span = [1 2];

Prompt(end+1,:) = {'Individual Bands','spectralbandsI'};
Formats(end+1,1).type = 'check';
Formats(end,1).format = 'input';
Formats(end,1).callback = {@do_spectralbandsI,'@ALL','@ALL'};
Formats(end,1).span = [1 2];

Prompt(end+1,:) = {'Filter','filter'};
Formats(end+1,1).type = 'list';
Formats(end,1).style = 'popupmenu';
Formats(end,1).format = 'input';
Formats(end,1).items = {'input','firls','butter','freq','wavelet','cheby','freq-shift'};
Formats(end,1).callback = {@correct_filter,'freqs','filter','freqs',header};
Formats(end,1).span = [1 3];

Formats(end+1,1).type = 'none';
Formats(end,1).span = [1 4];

Prompt(end+1,:) = {'Connectivity measures','measure'};
Formats(end+1,1).type = 'list';
Formats(end,1).style = 'listbox';
Formats(end,1).format = 'input';
Formats(end,1).items = {'PLI';'dPLI';'wPLI';'PLV';'wPLV';'PLT';'PTE';'rPTE';'SL';'SLc';'EC';'DTF';'TE'};
Formats(end,1).limits = [0 2];
Formats(end,1).size = [70 195];
Formats(end,1).callback = {@lab_get_connectivity,'@ALL','@ALL',data,header,cfg};
Formats(end,1).span = [5 2];

Prompt(end+1,:) = {'Marker','MARKER'};
Formats(end+1,1).type = 'check';
Formats(end,1).callback = {@lab_get_MARKERS,'@ALL','@ALL',header,cfg};

Prompt(end+1,:) = {'Detect Epochs','EPOCH'};
Formats(end+1,1).type = 'check';
Formats(end,1).callback = {@set_epochs,'EPOCH','EPOCH',header,cfg};

Prompt(end+1,:) = {'Phase','PHASE'};
Formats(end+1,1).type = 'check';
Formats(end,1).callback = {@lab_get_PHASE,'@ALL','@ALL'};

Prompt(end+1,:) = {'EC','EC'};
Formats(end+1,1).type = 'check';
Formats(end,1).callback = {@lab_get_EC,'@ALL','@ALL'};

Prompt(end+1,:) = {'SL','SL'};
Formats(end+1,1).type = 'check';
Formats(end,1).callback = {@lab_get_SL,'@ALL','@ALL'};

Prompt(end+1,:) = {'SLc','SLc'};
Formats(end+1,1).type = 'check';
Formats(end,1).callback = {@lab_get_SLc,'@ALL','@ALL'};

Prompt(end+1,:) = {'DTF','DTF'};
Formats(end+1,1).type = 'check';
Formats(end,1).callback = {@lab_get_DTF,'@ALL','@ALL',data,header};

Prompt(end+1,:) = {'TE','TE'};
Formats(end+1,1).type = 'check';
Formats(end,1).callback = {@lab_get_TE,'@ALL','@ALL'};

Formats(end+1,1).type = 'none';
Formats(end,1).span = [1 4];

Prompt(end+1,:) = {'Average matrix per patient','doaverage'};
Formats(end+1,1).type = 'check';
Formats(end,1).span = [1 4];

Formats(end+1,1).type = 'none';
Formats(end,1).span = [1 4];

Prompt(end+1,:) = {'Calculate phase randomization','randphase'};
Formats(end+1,1).type = 'check';
Formats(end,1).format = 'integer';
Formats(end,1).callback = {@set_randphase,'@ALL','@ALL'};
Formats(end,1).span = [1 2];

Prompt(end+1,:) = {'Iterations','numrandphase'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'integer';
Formats(end,1).limits = [0 9999];
Formats(end,1).size = 50;

Formats(end+1,1).type = 'none';
Formats(end,1).span = [1 4];

Prompt(end+1,:) = {'Correct for nearest neighboors (distance / 0=off)','correctdistance'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'float';
Formats(end,1).limits = [0 999];
Formats(end,1).size = 40;
Formats(end,1).span = [1 4];

Formats(end+1,1).type = 'none';
Formats(end,1).span = [1 4];

Prompt(end+1,:) = {'Write matrices','writematrix'};
Formats(end+1,1).type = 'check';
Formats(end,1).span = [1 4];

Prompt(end+1,:) = {'Plot matrices','PLOT'};
Formats(end+1,1).type = 'check';
Formats(end,1).callback = {@lab_set_plot_matrices,'PLOT','PLOT',header};
Formats(end,1).span = [1 4];

Prompt(end+1,:) = {'Graph analysis','GRAPH'};
Formats(end+1,1).type = 'check';
Formats(end,1).span = [1 3];
Formats(end,1).callback = {@set_graphanalysis,'@ALL','@ALL',cfg,header};

Prompt(end+1,:) = {'Kmeans clustering (maxclusters / 0=off)','clustering'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'integer';
Formats(end,1).limits = [0 999];
Formats(end,1).size = 50;
Formats(end,1).span = [1 3];

if nodel == 0
    Formats(end+1,1).type = 'none';
    Formats(end,1).span = [1 4];
    
    Prompt(end+1,:) = {'Delete results from previous run','deleteold'};
    Formats(end+1,1).type = 'check';
end

[cfg.CONNECT,Cancelled] = inputsdlg(Prompt,'Connectivity settings',Formats,cfg.CONNECT);
if isempty(cfg.CONNECT) | Cancelled == 1
    cfg.CONNECT = [];
    skipprocessing = 1;
    return
else
    pause(0.2);
end

if isfield(cfg.CONNECT,'GRAPH') & ~isempty(cfg.CONNECT.GRAPH) & isfield(cfg.CONNECT,'deleteold')
    cfg.GRAPH.deleteold = cfg.CONNECT.deleteold;
end

end

function freqs = correct_filter(filter,freqs,header)
   if strcmp(filter,'input')
       freqs = [0 0];
       if isfield(header,'highpass') & ~isempty(header.highpass)
           freqs(1,1) = header.highpass;
       end
       if isfield(header,'lowpass') & ~isempty(header.lowpass)
           freqs(1,2) = header.lowpass;
       end
   elseif isempty(freqs) | (size(freqs,1)== 1 & min(freqs) == 1 & max(freqs) == 48)
       freqs = [1 4;4 8;8 10;10 13;13 30;30 48];
   end
end

function settings = set_randphase(settings)
   if settings.randphase == false
       settings.randphase = true;
       if isempty(settings.numrandphase)
           settings.numrandphase = 100;
       end
   else
       settings.randphase = false;
       settings.numrandphase = [];
   end
end

function settings = set_graphanalysis(settings,cfg,header)
   if settings.doaverage == true
       if isfield(settings.GRAPH,'maxmatrix')
           settings.GRAPH.maxmatrix = 1;
       end
       if isfield(settings.GRAPH,'avgmatrix')
           settings.GRAPH.avgmatrix = false;
           settings.GRAPH.numavg = [];
       end
       settings = lab_set_graphanalysis(settings,1,cfg,header);
   else
       if isfield(settings.GRAPH,'maxmatrix')
           settings.GRAPH.maxmatrix = [];
       end
       settings = lab_set_graphanalysis(settings,[],cfg,header);
   end
end

function EPOCH = set_epochs(EPOCH,header,cfg)
    cfg.EPOCH = EPOCH;
    if isfield(cfg,'MICROST') & ~isempty(cfg.MICROST)
        cfg = lab_set_save_epochs(cfg,header,1,0,0,1);
    elseif isfield(header,'CORR') & ~isempty(header.CORR)
        cfg = lab_set_save_epochs(cfg,header,1,0,0,1);
    else
        cfg = lab_set_save_epochs(cfg,header,1,0,1,0);
    end
    EPOCH = cfg.EPOCH;
end

function settings = do_spectralbandsI(settings)
    if settings.spectralbandsI == false
        settings.spectralbandsI = true;
        settings.spectralbands = [2,4,5,6,8,9];
    else
        settings.spectralbandsI = false;
        settings.spectralbands = lab_get_spectralbands;
        settings.spectralbands = settings.spectralbands(2:end,:);
    end
end
