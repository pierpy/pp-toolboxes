% Read Cartool .rois or atlas(.hdr / eg AAL-atlas)
%
% [ROIS] = lab_create_rois(ROIS_file,ROISsettings,SPI,SPI_file)
%
% ROIS_file    = .rois file
%                .hdr atlas file (atlas information will be converted to rois
%                information using 'locs'
% ROISsettings = (optional)
% SPI          = structures variable with location info of solutionpoints
%                (=lab_read_spi); only needed for atlas conversion
% SPI_file     = Filename of MRI-file used for creation of solutionpoints
%
% written by F. Hatz 2012

function [ROIS,ROISsettings] = lab_create_rois(ROIS_file,ROISsettings,SPI,SPI_file)

if ~exist('ROISsettings','var')
    ROISsettings = [];
end

if exist('SPI','var') & ischar(SPI) & exist(SPI,'file')
    [~,~,~,SPI_format] = lab_filename(SPI);
    switch SPI_format
        case 'spi'
            SPI_file = [];
            SPI = lab_read_spi(SPI);
        case {'hdr','nii'}
            SPI_file = SPI;
            clearvars SPI
            SPI = lab_create_sp(SPI_file);
        otherwise
            ROIS = [];
            disp('    Abort, creating ROIS not possible, wrong SPI-input')
    end
elseif ~exist('SPI','var') | (~isfield(SPI,'x') & ~isfield(SPI,'anatomy'))
    [SPI_file,SPI_filepath]=uigetfile('*.spi;*.nii;*.hdr','Select SPI-file');
    [~,~,SPI_format] = lab_filename(SPI_file);
    if strcmp(SPI_format,'spi')
        SPI = lab_read_spi(fullfile(SPI_filepath,SPI_file));
    else
        SPI = lab_create_sp(fullfile(SPI_filepath,SPI_file));
    end
    SPI_file = fullfile(SPI_filepath,SPI_file);
    clearvars SPI_filepath
end

if ~exist('ROIS_file','var') | isempty(ROIS_file)
    [ROIS_file,ROIS_filepath]=uigetfile('*.rois;*.hdr;*.nii;*.els;*.xyz','Select ROIS-file');%G Bogaarts edit, allow electrode files
    if isempty(ROIS_file) | ROIS_file == 0
        ROIS = [];
        return
    else
        ROIS_file = fullfile(ROIS_filepath,ROIS_file);
    end
    clearvars ROIS_filepath
end

if isempty(ROIS_file) | ~exist(ROIS_file,'file')
    ROIS = [];
    return
end

if exist('SPI_file','var') & ~isempty(SPI_file) & ischar(SPI_file)
    [~,SPI_filepath,~,SPI_fileS] = lab_filename(SPI_file);
    ROISfileout = fullfile(SPI_filepath,[SPI_fileS '.rois']);
    clearvars SPI_filepath SPI_fileS
elseif ischar(ROIS_file)
    [~,TMP_filepath,~,TMP_fileS] = lab_filename(ROIS_file);
    ROISfileout = fullfile(TMP_filepath,[TMP_fileS '.rois']);
    clearvars TMP_filepath TMP_fileS
else
    ROISfileout = [];
end

if ischar(ROIS_file)
    [~,~,ROIS_format] = lab_filename(ROIS_file);
else
    ROIS_format = '';
end
if strcmp(ROIS_format,'rois')
    ROIS = lab_read_rois(ROIS_file);
    return
elseif strcmp(ROIS_format,'els') | strcmp(ROIS_format,'xyz') %G Bogaarts edit, Pseudo electrodes
    LOCS = lab_read_locs(ROIS_file);
    if (~exist('SPI','var') | isempty(SPI)) & exist('SPI_file','var') & exist(SPI_file,'file')
        SPI = lab_read_spi(SPI_file);
    end
    if exist('SPI','var') & ~isempty(SPI)
        ROIS = lab_create_rois_pseudo_electrodes(LOCS,SPI);
    else
        ROIS = [];
    end
    return
elseif exist(ROISfileout,'file')
    disp(['   Read ROIS from previous run (' lab_filename(ROISfileout) ')'])
    ROIS = lab_read_rois(ROISfileout);
    return
elseif exist('SPI','var') & ischar(ROIS_file) & exist(ROIS_file,'file')
    try
        atlas = lab_read_mri(ROIS_file);
    catch %#ok<CTCH>
        disp('    Creating ROIS not possible, ROIS-input should be atlas-file')
        ROIS = [];
        return
    end
elseif exist('SPI','var') & isstruct(ROIS_file) & isfield(ROIS_file,'anatomy')
    atlas = ROIS_file;
else
    disp('    Creating ROIS not possible, solutionpoints missing or no atlas-file for ROIS defined')
    ROIS = [];
    return
end

if ~isfield(ROISsettings,'CoregAtlas')
    button = questdlg('     coregister atlas to mri','Coregister atlas to mri','No','Yes','Yes');
    if strcmp(button,'Yes')
        ROISsettings.CoregAtlas = true;
    else
        ROISsettings.CoregAtlas = false;
    end
    button = questdlg('     remove border of regions','Remove border of regions','No','Yes','Yes');
    if strcmp(button,'Yes')
        ROISsettings.RemoveBorder = true;
    else
        ROISsettings.RemoveBorder = false;
    end
    pause(0.2);
end

if isfield(ROISsettings,'CoregAtlas') & ROISsettings.CoregAtlas == true
    if ~exist('SPI_file','var') | ~exist(SPI_file,'file')
        SPI_file = [];
    end
    if ~isfield(ROISsettings,'AtlasTemplate')
        ROISsettings.AtlasTemplate = '';
    end
    [atlas,atlasfile] = lab_match_atlas2mri(ROIS_file,SPI_file,ROISsettings.AtlasTemplate);
elseif ischar('ROIS_file') & exist(ROIS_file,'file')
    atlasfile = ROIS_file;
end
if isfield(ROISsettings,'RemoveBorder') & ROISsettings.RemoveBorder == true
    if exist('atlasfile','var') & ~isempty(atlasfile)
        [~,Atlas_filepath,~,Atlas_fileS] = lab_filename(atlasfile);
        if exist(fullfile(Atlas_filepath,[Atlas_fileS '_rmborder.hdr']),'file')
            atlas = lab_read_mri(fullfile(Atlas_filepath,[Atlas_fileS '_rmborder.hdr']));
        else
            atlas = lab_remove_borderzones(atlas);
            lab_write_hdr(fullfile(Atlas_filepath,[Atlas_fileS '_rmborder.hdr']),atlas);
        end
    else
        atlas = lab_remove_borderzones(atlas);
    end
end
roislist = unique(atlas.anatomy);
roislist = roislist(roislist > 0);
if length(roislist) == 116 & ~isfield(ROISsettings,'ReduceAAL')
    button = questdlg('Reduce AAL regions','Reduce regions','No','Yes','Yes');
    if strcmp(button,'Yes')
        ROISsettings.ReduceAAL = true;
    else
        ROISsettings.ReduceAAL = false;
    end
end
if length(roislist) == 116
    if ROISsettings.ReduceAAL == true
        [include,includetext] = lab_get_AALselection;
        disp(['     atlas is AAL, reduce to ' num2str(length(include)) ' regions (standard ' includetext ')'])
    else
        disp('     atlas is AAL')
        include = 1:116;
    end
    roislist = roislist(include);
    [tmp,tmp2] = lab_get_aal;
    ROIS.labels = tmp(include,1)';
    ROIS.label = tmp2(include,1)';
    clearvars tmp tmp2 include
else
    for i = 1:length(roislist)
        ROIS.labels{1,i} = ['ROI_' num2str(i)];
    end
end
ROIS.numrois = length(roislist) ;
ROIS.numsolutionpts = size(SPI.x,2);
for i = 1:size(SPI.x,2)
    x = round(SPI.x(i) + abs(atlas.transform(1,4)));
    y = round(SPI.y(i) + abs(atlas.transform(2,4)));
    z = round(SPI.z(i) + abs(atlas.transform(3,4)));
    locslabel(1,i) = atlas.anatomy(x,y,z);
end
clearvars x y z
ROIS.solutionptsAll = [];
for i =1:length(roislist)
    tmp = find(locslabel == roislist(i));
    ROIS.solutionpts{1,i} = tmp;
    ROIS.solutionptsAll = union(ROIS.solutionptsAll,tmp);
    clearvars tmp
end
if size(ROIS.solutionptsAll,1) > 1
    ROIS.solutionptsAll = ROIS.solutionptsAll';
end
ROIS.version = 'RO01';
if exist('ROISfileout','var')
    lab_write_rois(ROISfileout,ROIS);
end


% Create short label
ROIS = lab_create_rois_shortlabel(ROIS);