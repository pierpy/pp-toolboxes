function [cfg,skipprocessing] = lab_set_forwardsolution(cfg,header)

disp ('   Ask for forward solution-settings')

skipprocessing = 0;

if ~exist('cfg','var')
    cfg = [];
end

global THeader
if ~isempty(THeader)
    header = THeader;
elseif ~exist('header','var')
    header = [];
end

isfiles = lab_find_isfiles(cfg,header,[],true);

if ~isfield(cfg,'FWS') | ~isfield(cfg.FWS,'IS_file')
    cfg.FWS.folder = 'FWS';
    cfg.FWS.IS_file = isfiles.IS_file;
    cfg.FWS.SPI_file = isfiles.SPI_file;
    cfg.FWS.LOC_file = isfiles.LOC_file;
    cfg.FWS.LF = [];
    cfg.FWS.MRI = [];
    cfg.FWS.eformat = {'sef'};
    cfg.FWS.NOISE.mode = 3;
    cfg.FWS.NOISE.dB = 10;
    cfg.FWS.NOISE.coeff = 0;
    cfg.FWS.refchan = [];
end

Prompt = cell(0,2);
Formats = [];

Prompt(end+1,:) = {'Output-folder', 'folder'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'text';
Formats(end,1).size = 100;
Formats(end,1).span = [1 3];

Prompt(end+1,:) = {'IS-file','IS_file'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'file';
Formats(end,1).items = {'*.bin','Leadfield (.bin)';'*.mat','Headmodel (.mat)';'*.hdr','MRI-file (.hdr)'};
Formats(end,1).limits = [0 1]; % single file get
Formats(end,1).size = [400 0];

Prompt(end+1,:) = {'LF-Settings','LF'};
Formats(end+1,1).type = 'check';
Formats(end,1).callback = {@lab_get_LF,'LF','LF','IS_file','LOC_file'};

Prompt(end+1,:) = {'MRI-Settings','MRI'};
Formats(end+1,1).type = 'check';
Formats(end,1).callback = {@lab_get_pMRI,{'MRI','LF'},'MRI','IS_file','LF'};

Prompt(end+1,:) = {'Solutionpoints-file','SPI_file'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'file';
Formats(end,1).items = {'*.spi','Solutionpoints (.spi)';'*.hdr','MRI-file (.hdr)'};
Formats(end,1).limits = [0 1]; % single file get
Formats(end,1).size = [-1 0];

Prompt(end+1,:) = {'SPI-Settings','SPI'};
Formats(end+1,1).type = 'check';
Formats(end,1).callback = {@lab_get_spi,{'SPI','SPI_file'},'SPI','SPI_file','IS_file'};

Prompt(end+1,:) = {'ROIS','ROIS'};
Formats(end+1,1).type = 'check';
Formats(end,1).callback = {@lab_load_ROIS,'ROIS','ROIS','SPI_file','SPI'};

Prompt(end+1,:) = {'Electrodes-file','LOC_file'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'file';
Formats(end,1).items = {'*.els;*.xyz','Electrodes-file (.els/.xyz)'};
Formats(end,1).limits = [0 1]; % single file get
Formats(end,1).size = [-1 0];

Prompt{end+1,1} = ' ';
Formats(end+1,1).type = 'text';
Formats(end,1).span = [1 2];

Prompt(end+1,:) = {'Reference channel(s)','refchan'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'vector';
Formats(end,1).limits = [-inf inf];
Formats(end,1).size = 80;
Formats(end,1).span = [1 3];

Prompt(end+1,:) = {'Add noise','NOISE'};
Formats(end+1,1).type = 'check';
Formats(end,1).callback = {@lab_set_add_noise,'@ALL','@ALL'};
Formats(end,1).span = [1 3];

Prompt(end+1,:) = {'File format','eformat'};
Formats(end+1,1).type = 'list';
Formats(end,1).style = 'listbox';
Formats(end,1).format = 'input';
Formats(end,1).items = {'sef';'edf';'eph';'ep';'txt'};
Formats(end,1).limits = [0 2]; % multi-select
Formats(end,1).size = [60 85];
Formats(end,1).callback = {@lab_get_format,'scaletxt','format','scaletxt'};

[cfg.FWS,Cancelled] = inputsdlg(Prompt,'Forward solution',Formats,cfg.FWS);
if isempty(cfg.FWS) | Cancelled == 1
    cfg.FWS = [];
    skipprocessing = 1;
    return
end

% read *.spi file
if ~isfield(cfg.FWS,'locs') & exist(cfg.FWS.SPI_file,'file')
    disp('   Read SPI')
    if strcmp(cfg.FWS.SPI_file(end-3:end),'.spi')
        locs = lab_read_spi(cfg.FWS.SPI_file);
    else
        [locs,cfg.FWS] = lab_create_sp(cfg.FWS.SPI_file,cfg.FWS);
    end
    cfg.FWS.locs = locs;
    cfg.FWS.numsp = size(cfg.FWS.locs.x,2);
end
% read electrodes file
if ~isfield(cfg.FWS,'elec') & exist(cfg.FWS.LOC_file,'file')
    disp('   Read LOCS')
    cfg.FWS.elec = lab_read_locs(cfg.FWS.LOC_file);
    cfg.FWS.numelec = size(cfg.FWS.elec.x,2);
end
% read LF-bin
if ~isfield(cfg.FWS,'LFbin') & exist(cfg.FWS.IS_file,'file') & strcmp(cfg.FWS.IS_file(end-3:end),'.bin') & isfield(cfg.FWS,'numelec')  & isfield(cfg.FWS,'numsp')
    disp('   Read LFbin')
    cfg.FWS.LFbin = lab_read_LFbin(cfg.FWS.IS_file,cfg.FWS.numelec,cfg.FWS.numsp);
end

end