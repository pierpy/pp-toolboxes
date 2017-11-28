function [cfg,skipprocessing] = lab_set_dipolfit(cfg,header,data,doepochs)
disp ('   Ask for inverse-solution-settings')
skipprocessing = 0;

if ~exist('doepochs','var')
    doepochs = 0;
end
if ~exist('data','var')
    data = [];
end
if ~exist('header','var')
    header = [];
end
if ~exist('cfg','var')
    cfg = [];
end

if isfield(cfg,'BADELEC') & isfield(cfg.BADELEC,'markbad') & cfg.BADELEC.markbad == true
    markerbad = 1; %#ok<NASGU>
else
    markerbad = 0; %#ok<NASGU>
end

% Search for IS files
isfiles = lab_find_isfiles(cfg,header);

if ~isfield(cfg,'IS') | ~isfield(cfg.IS,'IS_file')
    cfg.IS.averageepochs = false;
    cfg.IS.IS_file = ' ';
    cfg.IS.SPI_file = 'None';
    cfg.IS.SPI =[];
    cfg.IS.LOC_file = ' ';
    cfg.IS.MRI_file = ' ';
    cfg.IS.ROIS_file = ' ';
    cfg.IS.ROISmethods = {};
    cfg.IS.spectralbands = [];
    cfg.IS.spectralbandsI = false;
    cfg.IS.ROISsettings = [];
    cfg.IS.LF = [];
    cfg.IS.eegsource = 'input';
    cfg.IS.LAPL = [];
    cfg.IS.interpolatebad = true;
    cfg.IS.eformat = {'ris'};
    cfg.IS.resultscalar = true;
    cfg.IS.IndividualFiles = false;
    cfg.IS.type = 'Dipol fit';
    cfg.IS.issettings.numdipols = 1;
    cfg.IS.issettings.model = 'regional';
    cfg.IS.issettings.symmetry = '';
end
if strcmp(cfg.IS.SPI_file,' ') | isempty(cfg.IS.SPI_file)
    cfg.IS.SPI_file = 'None';
end

cfg.IS = lab_ISconvertnames(cfg.IS);

Prompt = cell(0,2);
Formats = [];

if doepochs == 1
    Prompt(end+1,:) = {'Average epochs','averageepochs'};
    Formats(end+1,1).type = 'check';
    Formats(end,1).span = [1 3];
    
    Formats(end+1,1).type = 'none';
    Formats(end,1).span = [1 3];
end

Prompt(end+1,:) = {'Reference','eegsource'};
Formats(end+1,1).type = 'list';
Formats(end,1).style = 'popupmenu';
Formats(end,1).format = 'input';
Formats(end,1).items = {'input';'mean';'median';'laplacian'};
Formats(end,1).callback = {@lab_get_eegsource,'@ALL','@ALL',cfg,header};

Prompt(end+1,:) = {'Laplacian','LAPL'};
Formats(end+1,1).type = 'check';
Formats(end,1).callback = {@lab_get_laplacian,'LAPL','LAPL'};
Formats(end,1).span = [1 2];

Formats(end+1,1).type = 'none';
Formats(end,1).span = [1 3];

Prompt(end+1,:) = {'Input-file for Inverse solution','IS_file'};
Formats(end+1,1).type = 'list';
Formats(end,1).style = 'popupmenu';
Formats(end,1).format = 'input';
Formats(end,1).items = {'Select File','Leadfield (LF.bin)', ...
    'Headmodel (.mat)','MRI-file (.hdr)','MRI-file (.nii)','MRI-file (iso.fif)', ...
    'MRI-file (mri.fif)','MRI-file (dicom)'};
Formats(end,1).size = 300;
Formats(end,1).callback = {@lab_get_isfile,'@ALL','@ALL',isfiles,0};
Formats(end,1).span = [1 3];

Prompt(end+1,:) = {'Individual Files','IndividualFiles'};
Formats(end+1,1).type = 'check';
Formats(end,1).enable = 'inactive';

Prompt(end+1,:) = {'MRI settings','MRI'};
Formats(end+1,1).type = 'check';
Formats(end,1).callback = {@lab_get_pMRI,{'MRI','LF'},'MRI','IS_file','LF',isfiles};

Prompt(end+1,:) = {'Leadfield','LF'};
Formats(end+1,1).type = 'check';
Formats(end,1).callback = {@lab_get_LF,'LF','LF','IS_file','LOC_file'};

Formats(end+1,1).type = 'none';
Formats(end,1).span = [1 3];

Prompt(end+1,:) = {'Grid search','SPI_file'};
Formats(end+1,1).type = 'list';
Formats(end,1).style = 'popupmenu';
Formats(end,1).format = 'input';
Formats(end,1).items = {'None','Select File','Solutionpoints (.spi)','MRI-file (IS-file)','MRI-file (spi.hdr)','MRI-file (spi.nii)'};
Formats(end,1).size = 250;
Formats(end,1).callback = {@lab_get_spi,{'SPI','SPI_file'},'SPI','SPI_file','IS_file'};
Formats(end,1).span = [1 2];

Prompt(end+1,:) = {'Create grid','SPI'};
Formats(end+1,1).type = 'check';
Formats(end,1).callback = {@lab_get_spi,{'SPI','SPI_file'},'SPI','SPI_file','IS_file'};
% Formats(end,1).span = [1 3];

Formats(end+1,1).type = 'none';
Formats(end,1).span = [1 3];

Prompt(end+1,:) = {'Type of Electrodes file','LOC_file'};
Formats(end+1,1).type = 'list';
Formats(end,1).style = 'popupmenu';
Formats(end,1).format = 'input';
Formats(end,1).items = {'Select File','Locations in input file','Electrodes-file (.xyz)','Electrodes-file (.els)'};
Formats(end,1).size = 300;
Formats(end,1).callback = {@lab_get_LOC,'LOC_file','LOC_file'};
Formats(end,1).span = [1 3];

Formats(end+1,1).type = 'none';
Formats(end,1).span = [1 3];

Prompt(end+1,:) = {'MRI-file (visualisation)','MRI_file'};
Formats(end+1,1).type = 'list';
Formats(end,1).style = 'popupmenu';
Formats(end,1).format = 'input';
Formats(end,1).items = {'Select File','MRI-file (IS-file)','MRI-file (.hdr)'};
Formats(end,1).size = 300;
Formats(end,1).callback = {@lab_get_MRI,'MRI_file','MRI_file'};
Formats(end,1).span = [1 3];

Formats(end+1,1).type = 'none';
Formats(end,1).span = [1 3];
    
Prompt(end+1,:) = {'File format','eformat'};
Formats(end+1,1).type = 'list';
Formats(end,1).style = 'listbox';
Formats(end,1).format = 'input';
Formats(end,1).items = {'ris';'sef';'edf';'eph';'ep';'txt'};
Formats(end,1).limits = [0 2]; % multi-select
Formats(end,1).size = [60 90];
Formats(end,1).callback = {@lab_get_format,'scaletxt','eformat','scaletxt'};
Formats(end,1).span = [2 1];

Prompt(end+1,:) = {'Scalar result','resultscalar'};
Formats(end+1,1).type = 'check';
Formats(end,1).span = [1 2];

Prompt(end+1,:) = {'Scale TXT','scaletxt'};
Formats(end+1,1).type = 'check';
Formats(end,1).callback = {@lab_get_format,'scaletxt','eformat','scaletxt'};
Formats(end,1).span = [1 2];

[cfg.IS,Cancelled] = inputsdlg(Prompt,'Dipol fit',Formats,cfg.IS);
if isempty(cfg.IS) | Cancelled == 1
    cfg.IS = [];
    skipprocessing = 1;
    return
else
    cfg.IS = lab_ISconvertnames2(cfg.IS);
    if strcmp(cfg.IS.SPI_file,'None')
        cfg.IS.SPI_file = '';
    end
end

% read SPI-file
if ~isfield(cfg.IS,'locs') & exist(cfg.IS.SPI_file,'file')
    disp('   Read SPI')
    if strcmp(cfg.IS.SPI_file(end-3:end),'.spi')
        locs = lab_read_spi(cfg.IS.SPI_file);
    else
        [locs,cfg.IS] = lab_create_sp(cfg.IS.SPI_file,cfg.IS.SPI);
    end
    cfg.IS.locs = locs;
end
% read IS-file
if ~isempty(data)
    cfg.IS = lab_read_ismatrix(cfg.IS,cfg,data);
else
    cfg.IS = lab_read_ismatrix(cfg.IS,cfg);
end
cfg.IS.COV = [];
