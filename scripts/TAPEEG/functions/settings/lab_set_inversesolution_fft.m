function [cfg,skipprocessing] = lab_set_inversesolution_fft(cfg,header,data)

disp ('   Ask for inverse-solution-fft-settings')

skipprocessing = 0;

if ~exist('cfg','var')
    cfg = [];
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

% Search for IS files
isfiles = lab_find_isfiles(cfg,header);

if ~isfield(cfg,'ISFFT') | ~isfield(cfg.ISFFT,'IS_file')
    cfg.ISFFT.folder = 'SpectIS';
    cfg.ISFFT.IS_file = ' ';
    cfg.ISFFT.SPI_file = ' ';
    cfg.ISFFT.SPI =[];
    cfg.ISFFT.LOC_file = ' ';
    cfg.ISFFT.MRI_file = ' ';
    cfg.ISFFT.ROIS_file = ' ';
    cfg.ISFFT.ROISmethods = {};
    cfg.ISFFT.ROISsettings = [];
    cfg.ISFFT.LF = [];
    cfg.ISFFT.COV = [];
    cfg.ISFFT.IndividualFiles = false;
    cfg.ISFFT.eegsource = 'input';
    cfg.ISFFT.interpolatebad = true;
    cfg.ISFFT.excludebad = false;
    cfg.ISFFT.eformat = {'ris'};
    cfg.ISFFT.resultscalar = true;
    cfg.ISFFT.deleteold = true;
    cfg.ISFFT.FFT.winsize = 4;
    cfg.ISFFT.FFT.hanpercent = 20;
    cfg.ISFFT.FFT.lowfreq = 1;
    cfg.ISFFT.FFT.highfreq = 50;
    cfg.ISFFT.FFT.spectralbands = lab_get_spectralbands;
    cfg.ISFFT.EPOCH.minimalpart = 0.4;
    cfg.ISFFT.EPOCH.percentgood = 100;
    cfg.ISFFT.EPOCH.markerexclude = {''};
    cfg.ISFFT.EPOCH.markerinclude = {''};
    cfg.ISFFT.EPOCH.markerstart ='';
    cfg.ISFFT.EPOCH.markerstop = '';
    cfg.ISFFT.EPOCH.interpolate3D = true;
    cfg.ISFFT.EPOCH.BAD.freqlim50 = 50;
    cfg.ISFFT.EPOCH.BAD.freqlim60 = 50;
    cfg.ISFFT.EPOCH.BAD.freqlimlow = 70;
    cfg.ISFFT.EPOCH.BAD.freqlimhigh = 50;
    cfg.ISFFT.EPOCH.BAD.spectslow = [0 2];
    cfg.ISFFT.EPOCH.BAD.spectshigh = [15 50];
    cfg.ISFFT.EPOCH.BAD.zvaluebroken = 0;
    cfg.ISFFT.EPOCH.BAD.zvaluevars = 3;
    cfg.ISFFT.EPOCH.BAD.zvaluehurst = 3;
    cfg.ISFFT.EPOCH.BAD.zvaluemedian = 3;
    cfg.ISFFT.EPOCH.BAD.zvaluecorr = 3;
    cfg.ISFFT.EPOCH.BAD.LAPL.lap_maxdistance = 2.5;
    cfg.ISFFT.EPOCH.BAD.LAPL.lap_weightmaxdistance = 30;
    cfg.ISFFT.EPOCH.BAD.LAPL.lap_excluderef = true;
    cfg.ISFFT.EPOCH.BAD.PEAK2MIN.lowfreqpeak = 4;
    cfg.ISFFT.EPOCH.BAD.PEAK2MIN.highfreqpeak = 14;
    cfg.ISFFT.EPOCH.BAD.PEAK2MIN.MinPeak2Min = 1.3;
    cfg.ISFFT.EPOCH.BAD.PEAK2MIN.mode = 'weighted';
    cfg.ISFFT.EPOCH.BAD.PEAK2MIN.threshold = [];
    cfg.ISFFT.EPOCH.BAD.PEAK2MIN.factor = 3;
    if exist('header','var') & isfield(header,'events') & isfield(header.events,'TYP') & ~isempty(header.events.TYP)
        cfg.ISFFT.EPOCH.TRIGGER = unique(header.events.TYP);
    end
end
if isfield(cfg.ISFFT.FFT,'spectralbands') & isnumeric(cfg.ISFFT.FFT.spectralbands) & size(cfg.ISFFT.FFT.spectralbands,2) == 2
    cfg.ISFFT.FFT.spectralbands = cat(2,cell(size(cfg.ISFFT.FFT.spectralbands,1),1),num2cell(cfg.ISFFT.FFT.spectralbands));
end

cfg.ISFFT = lab_ISconvertnames(cfg.ISFFT);

Prompt = cell(0,2);
Formats = [];

Prompt(end+1,:) = {'Output-folder', 'folder'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'text';
Formats(end,1).size = 100;
Formats(end,1).span = [1 3];

Prompt(end+1,:) = {'Reference','eegsource'};
Formats(end+1,1).type = 'list';
Formats(end,1).style = 'popupmenu';
Formats(end,1).format = 'input';
Formats(end,1).items = {'input';'mean';'median';'laplacian'};
Formats(end,1).callback = {@lab_get_eegsource,'@ALL','@ALL',cfg,header};

Prompt(end+1,:) = {'Laplacian','LAPL'};
Formats(end+1,1).type = 'check';
Formats(end,1).callback = {@lab_get_laplacian,'LAPL','LAPL'};

Prompt(end+1,:) = {'Interpolate bad channels','interpolatebad'};
Formats(end+1,1).type = 'check';

Formats(end+1,1).type = 'none';
Formats(end,1).span = [1 2];

Prompt(end+1,:) = {'Exclude bad channels','excludebad'};
Formats(end+1,1).type = 'check';
Formats(end,1).format = 'integer';
Formats(end,1).enable = 'inactive';
Formats(end,1).callback = {@set_excludebad,'@ALL','@ALL'};

Prompt(end+1,:) = {'EPOCH-settings','EPOCH'};
Formats(end+1,1).type = 'check';
Formats(end,1).callback = {@set_epochs,'EPOCH','EPOCH',header,cfg};
Formats(end,1).span = [1 3];

Prompt(end+1,:) = {'FFT-settings','FFT'};
Formats(end+1,1).type = 'check';
Formats(end,1).callback = {@lab_set_calculate_spectras,'FFT','FFT',header,1};
Formats(end,1).span = [1 3];

Formats(end+1,1).type = 'none';
Formats(end,1).span = [1 3];

[cfg.ISFFT,Prompt,Formats] = lab_get_isfiles(cfg.ISFFT,0,1,isfiles,Prompt,Formats);

Prompt(end+1,:) = {'File format','eformat'};
Formats(end+1,1).type = 'list';
Formats(end,1).style = 'listbox';
Formats(end,1).format = 'input';
Formats(end,1).items = {'ris';'sef';'edf';'eph';'ep';'txt'};
Formats(end,1).limits = [0 2]; % multi-select
Formats(end,1).size = [60 90];
Formats(end,1).callback = {@lab_get_format,'scaletxt','eformat','scaletxt'};

Prompt(end+1,:) = {'Scalar result','resultscalar'};
Formats(end+1,1).type = 'check';
Formats(end,1).size = [-1 -1];

Prompt(end+1,:) = {'Scale TXT','scaletxt'};
Formats(end+1,1).type = 'check';
Formats(end,1).callback = {@lab_get_format,'scaletxt','eformat','scaletxt'};

Prompt(end+1,:) = {'Delete results from previous run','deleteold'};
Formats(end+1,1).type = 'check';
Formats(end,1).size = [-1 -1];

[cfg.ISFFT,Cancelled] = inputsdlg(Prompt,'Inverse solution FFT',Formats,cfg.ISFFT);
if isempty(cfg.ISFFT) | Cancelled == 1
    cfg.ISFFT = [];
    skipprocessing = 1;
    return
else
    cfg.ISFFT = lab_ISconvertnames2(cfg.ISFFT);
end

% read *.spi file
if ~isfield(cfg.ISFFT,'slocs') & exist(cfg.ISFFT.SPI_file,'file')
    disp('   Read SPI')
    if strcmp(cfg.ISFFT.SPI_file(end-3:end),'.spi')
        slocs = lab_read_spi(cfg.ISFFT.SPI_file);
    else
        [slocs,cfg.ISFFT] = lab_create_sp(cfg.ISFFT.SPI_file,cfg.ISFFT.SPI);
    end
    cfg.ISFFT.slocs = slocs;
end
% read rois
if ~isfield(cfg.ISFFT,'rois') & exist(cfg.ISFFT.ROIS_file,'file')
    disp('   Read Rois')
    [cfg.ISFFT.rois,cfg.ISFFT] = lab_create_rois(cfg.ISFFT.ROIS_file,cfg.ISFFT,cfg.ISFFT.slocs,cfg.ISFFT.SPI_file);
end
% read IS matrix
if exist('data','var')
    cfg.ISFFT = lab_read_ismatrix(cfg.ISFFT,cfg,data);
else
    cfg.ISFFT = lab_read_ismatrix(cfg.ISFFT,cfg);
end

if isfield(cfg.ISFFT,'LF') & (~isfield(cfg.ISFFT.LF,'locs') | isempty(cfg.ISFFT.LF.locs)) & isfield(cfg.ISFFT,'LOC_file') & exist(cfg.ISFFT.LOC_file,'file')
    cfg.ISFFT.LF.locs = lab_read_locs(cfg.ISFFT.LOC_file);
end

end

function settings = set_excludebad(settings)
   settings2 = lab_ISconvertnames2(settings);
   [~,~,ISformat] = lab_filename(settings2.IS_file);
   if strcmp(ISformat,'is') | strcmp(ISformat,'spinv')
       settings.excludebad = false;
       return
   end
   if settings.excludebad == false
       settings.excludebad = true;
       settings.eegsource = 'input';
   else
       settings.excludebad = false;
   end

end

function EPOCH = set_epochs(EPOCH,header,cfg)
    cfg.EPOCH = EPOCH;
    cfg = lab_set_save_epochs(cfg,header,1,1,1);
    EPOCH = cfg.EPOCH;
end

