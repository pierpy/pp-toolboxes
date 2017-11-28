% Create solutionpoints based on mri-information
%   Solutionpoints are placed on a grid in gray or gray&white matter /
%   surface of brain
%
%   [locs,settings] = lab_create_sp(mrifile,settings)
%         settings.method = 'Gray' / 'GrayWhite' / 'Surface'
%         settings.numspi = Minimal number of resulting solutionpoints
%         settings.atlasfile (optional) = mask mri-file by atlas
%         settings.atlasexclude = atlas regions to exclude
%
% written by F. Hatz 2012

function [locs,settings] = lab_create_sp(mrifile,settings)

if exist('settings','var') & ~isstruct(settings) & isnumeric(settings)
    tmp = settings;
    clearvars settings
    settings.numspi = tmp;
    clearvars tmp
end
if ~exist('settings','var')
    settings = [];
end
if ~isfield(settings,'SPI') | ~isfield(settings,'method')
    settings.method = 'Gray';
end
if ~exist('mrifile','var') | isempty(mrifile)
    [filename,filepath] = uigetfile('*.hdr;*.fif;*.nii','select MRI-file');
    mrifile = fullfile(filepath,filename);
    clearvars filename filepath
end

if ischar(mrifile)
    try
        [~,filename] = lab_read_mri(mrifile);
    catch %#ok<CTCH>
        locs = [];
        return
    end
    if isempty(filename) | ~exist(filename,'file')
        return
    end
    mri = load_nii(filename);
    
    [~,filepath,~,filenameS] = lab_filename(filename);
    warning off %#ok<WNOFF>
    mkdir(fullfile(filepath,'MriSeg'));
    warning on %#ok<WNON>
    filenameout2 = fullfile(fullfile(filepath,'MriSeg'),filenameS);
    filenameout = fullfile(filepath,filenameS);
    clearvars filepath filenameS
    if exist([filenameout '_locs.spi'],'file')
        disp('    read source locs from previous run')
        locs = lab_read_spi([filenameout '_locs.spi']);
        return
    end
    if ~strcmp(settings.method,'All')
        if ~exist([filenameout2 '_gray.hdr'],'file')
            settings.mrifile = filename;
            settings.SEGbrain = 1;
            settings.SEGprobmaps = 0;
            settings.SEGcorrect = 1;
            lab_segment_mri(settings);
        end
        mri = load_nii([filenameout2 '_gray.hdr']);
        if strcmp(settings.method,'GrayWhite') | strcmp(settings.method,'Surface')
            mri2 = load_nii([filenameout2 '_white.hdr']);
            mri.img(mri2.img == 1) = 1;
            clearvars mri2
        end
    end
else
    mri = mrifile;
    filenameout = [];
end
if isstruct(mri) & isfield(mri,'img')
    vol = mri.img;
elseif isstruct(mri) & isfield(mri,'anatomy')
    vol = mri.anatomy;
else
    locs = [];
    return
end

% Get originator coordinates
if isfield(mri,'hdr') & isfield(mri.hdr,'hist') & isfield(mri.hdr.hist,'originator')
    originator = mri.hdr.hist.originator(1,1:3);
    locs.orig = [0 0 0];
elseif isfield(mri,'transform')
    originator = abs(mri.transform(1:3,4)');
    locs.orig = [0 0 0];
else
    originator = [1 1 1];
end

if isfield(settings,'numspi') & settings.numspi > 0
    numspi = round(settings.numspi);
else
    numspi = 10000;
end

if isfield(settings,'atlasfile') & exist(settings.atlasfile,'file')
    if ~isfield(settings,'AtlasTemplate');
        settings.AtlasTemplate = '';
    end
    [atlas,atlasfile] = lab_match_atlas2mri(settings.atlasfile,mrifile,settings.AtlasTemplate);
    if isfield(settings,'RemoveBorder') & settings.RemoveBorder == true
        atlas = lab_remove_borderzones(atlas);
        if ~isempty(atlasfile)
            lab_write_hdr(atlasfile,atlas);
        end
    end
    exclude = [];
    roislist = unique(atlas.anatomy);
    roislist = roislist(roislist > 0);
    if isfield(settings,'atlasexclude') & ~isempty(settings.atlasexclude)
        disp(['     atlas: excluded regions' num2str(settings.atlasexclude)])
        exclude = roislist(settings.atlasexclude);
    elseif length(roislist) == 116 & isfield(settings,'ReduceAAL') & settings.ReduceAAL == true
        [includechannels,includetext] = lab_get_AALselection;
        disp(['     atlas is AAL, reduce to standard regions (standard ' includetext ')'])
        exclude = roislist(setdiff(1:length(roislist),includechannels));
    end
    if ~isempty(exclude)
        for i = 1:length(exclude)
            atlas.anatomy(atlas.anatomy==exclude(i)) = 0;
        end
    end
    vol(atlas.anatomy==0) = 0;
end

if ~strcmp(settings.method,'Surface')
    numsolutions = 0;
    scale = 10;
    while numsolutions < numspi
        % Scale down image
        xorig = round((size(vol,1)/scale-floor(size(vol,1)/scale))*scale/2);
        yorig = round((size(vol,2)/scale-floor(size(vol,2)/scale))*scale/2);
        zorig = round((size(vol,3)/scale-floor(size(vol,3)/scale))*scale/2);
        vol2 = vol(xorig+1:scale:end,yorig+1:scale:end,zorig+1:scale:end);
        % Get solutionpoints
        locs.x = [];
        locs.y = [];
        locs.z = [];
        for i = 1:size(vol2,3)
            [x,y] = find(vol2(:,:,i) > 0);
            if ~isempty(x)
                locs.x = [locs.x ((x-1)'*scale - originator(1,1) + xorig)];
                locs.y = [locs.y ((y-1)'*scale - originator(1,2) + yorig)];
                locs.z = [locs.z ((ones(1,size(x,1))*i-1)*scale - originator(1,3) + zorig)];
            end
        end
        locs.scale = scale;
        scale = scale - 1;
        numsolutions = size(locs.x,2);
    end
else
    cfg.nvert = numspi;
    mri.brain = vol;
    bnd = lab_mesh_segm(mri,cfg);
    locs.x = bnd.pnt(:,1)';
    locs.y = bnd.pnt(:,2)';
    locs.z = bnd.pnt(:,3)';
end

for i = 1:size(locs.x,2)
    locs.labels{1,i} = ['SP_' num2str(i)];
end
locs  = lab_locs2sph(locs);
if ~isempty(filenameout)
    lab_write_spi([filenameout '_locs.spi'],locs);
end