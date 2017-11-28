function [cfg,skipprocessing] = lab_set_inversesolution(cfg,header)
disp ('   Ask for inverse-solution-settings')

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

% Search for IS files
isfiles = lab_find_isfiles(cfg,header);

if ~isfield(cfg,'IS') | ~isfield(cfg.IS,'IS_file')
    cfg.IS.folder = 'IS';
    cfg.IS.IS_file = ' ';
    cfg.IS.SPI_file = ' ';
    cfg.IS.SPI =[];
    cfg.IS.LOC_file = ' ';
    cfg.IS.MRI_file = ' ';
    cfg.IS.ROIS_file = ' ';
    cfg.IS.ROISmethods = {};
    cfg.IS.spectralbands = lab_get_spectralbands;
    cfg.IS.spectralbandsI = false;
    cfg.IS.ROISsettings = [];
    cfg.IS.LF = [];
    cfg.IS.MRI = [];
    cfg.IS.COV.covmethod = 'Input-File';
    cfg.IS.eegsource = 'input';
    cfg.IS.interpolatebad = true;
    cfg.IS.excludebad = false;
    cfg.IS.eformat = {'ris'};
    cfg.IS.resultscalar = true;
    cfg.IS.IndividualFiles = false;
end

cfg.IS = lab_ISconvertnames(cfg.IS);

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

Formats(end+1,1).type = 'none';
Formats(end,1).span = [1 3];

[cfg.IS,Prompt,Formats] = lab_get_isfiles(cfg.IS,0,0,isfiles,Prompt,Formats);

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

Prompt(end+1,:) = {'Scale TXT','scaletxt'};
Formats(end+1,1).type = 'check';
Formats(end,1).callback = {@lab_get_format,'scaletxt','eformat','scaletxt'};

[cfg.IS,Cancelled] = inputsdlg(Prompt,'Inverse solution',Formats,cfg.IS);
if isempty(cfg.IS) | Cancelled == 1
    cfg.IS = [];
    skipprocessing = 1;
    return
else
    pause(0.2);
    cfg.IS = lab_ISconvertnames2(cfg.IS);
end

% read SPI-file
if isfield(cfg.IS,'SPI_file') & exist(cfg.IS.SPI_file,'file')
    disp('   Read SPI')
    if strcmp(cfg.IS.SPI_file(end-3:end),'.spi')
        slocs = lab_read_spi(cfg.IS.SPI_file);
    else
        [slocs,cfg.IS.SPI] = lab_create_sp(cfg.IS.SPI_file,cfg.IS.SPI);
    end
    cfg.IS.slocs = slocs;
end
% read ROIS-file
if isfield(cfg.IS,'ROIS_file') & exist(cfg.IS.ROIS_file,'file') & isfield(cfg.IS,'slocs')
    disp('   Read Rois')
    [rois,cfg.IS] = lab_create_rois(cfg.IS.ROIS_file,cfg.IS,cfg.IS.slocs,cfg.IS.SPI_file);
    cfg.IS.rois = rois;
end
% read IS-file
if exist('header','var') & isfield(header,'numdatachannels')
    cfg.IS = lab_read_ismatrix(cfg.IS,cfg,header.numdatachannels);
else
    cfg.IS = lab_read_ismatrix(cfg.IS,cfg);
end
% read LOC-file
if isfield(cfg.IS,'LF') & isfield(cfg.IS,'LOC_file') & exist(cfg.IS.LOC_file,'file')
    cfg.IS.LF.locs = lab_read_locs(cfg.IS.LOC_file);
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
