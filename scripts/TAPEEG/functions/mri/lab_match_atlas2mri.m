% Match atlas-template to mri-file with spm8 routines
%
% atlasout = lab_match_atlas2mri(atlas,mri)
%
%    atlas =    filename of atlas file or mri-structure (from lab_read_mri)
%    mri   =    filename of atlas file or mri-structure (from lab_read_mri)
%
%    atlasout = structure with coregistered atlas
%
% F. Hatz 2013

function [Atlasout,AtlasFileout] = lab_match_atlas2mri(atlas,mri,atlastemplate)

disp('   Match Atlas to MRI')
rng('default');
rng('shuffle');

if ~exist('atlas','var')
    [settings,skipprocessing] = lab_set_match_atlas2mri;
    if skipprocessing == 1
        Atlasout = [];
        AtlasFileout = [];
        return
    end
    atlas = settings.Atlas_file;
    mri = settings.MRI_file;
    atlastemplate = settings.AtlasTemplate_file;
elseif ~exist('atlastemplate','var')
    atlastemplate = '';
end

if exist('mri','var') & isstruct(mri) & isfield(mri,'anatomy')
    MriFileout = fullfile(pwd,'MRI.nii');
    ft_write_mri(MriFileout,int8(mri.anatomy),'transform',mri.transform,'dataformat','nifti');
else
    if ~exist('mri','var') | ~ischar(mri) | ~exist(mri,'file') | strcmp(mri(end-2:end),'spi')
        [tmp1,tmp2] = uigetfile('*.nii;*.hdr','Select MRI file');
        if tmp1 == 0
            Atlasout = [];
            AtlasFileout = [];
            return
        end
        mri = fullfile(tmp2,tmp1);
        clearvars tmp1 tmp2
    end
    [mri,mri_file] = lab_read_mri(mri);
    [mri_file,mri_filepath,~,mri_fileS] = lab_filename(mri_file);
    warning off %#ok<WNOFF>
    mkdir(fullfile(mri_filepath,'MriSeg'));
    warning on %#ok<WNON>
    MriFileout = fullfile(fullfile(mri_filepath,'MriSeg'),mri_fileS);
end

if exist('atlas','var') & isstruct(atlas) & isfield(atlas,'anatomy')
    if exist('mri_filepath','var')
        atlas_filepathout = fullfile(mri_filepath,'Coreg');
    else
        atlas_filepathout = pwd;
    end
    atlas_fileS = 'Atlas';
    AtlasFileout = fullfile(atlas_filepathout,atlas_fileS);
    ft_write_mri([AtlasFileout '.nii'],int16(atlas.anatomy),'transform',atlas.transform,'dataformat','nifti');
else
    if ~exist('atlas','var') | ~ischar(atlas) | ~exist(atlas,'file')
        [tmp1,tmp2] = uigetfile('*.nii;*.hdr','Select Atlas file');
        if tmp1 == 0
            Atlasout = [];
            AtlasFileout = [];
            return
        end
        atlas = fullfile(tmp2,tmp1);
        clearvars tmp1 tmp2
    end
    [atlas,atlas_file] = lab_read_mri(atlas);
    [atlas_file,atlas_filepath,~,atlas_fileS] = lab_filename(atlas_file);
    if exist('mri_filepath','var')
        atlas_filepathout = fullfile(mri_filepath,'Coreg');
    else
        atlas_filepathout = fullfile(atlas_filepath,'Coreg');
    end
    warning off %#ok<WNOFF>
    mkdir(atlas_filepathout);
    warning on %#ok<WNON>
    AtlasFileout = fullfile(atlas_filepathout,atlas_fileS);
    if ~exist([AtlasFileout '.nii'],'file')
        ft_write_mri([AtlasFileout '.nii'],int16(atlas.anatomy),'transform',atlas.transform,'dataformat','nifti');
        if exist(fullfile(atlas_filepath,[atlas_fileS '.txt']),'file')
            copyfile(fullfile(atlas_filepath,[atlas_fileS '.txt']),[AtlasFileout '.txt']);
        end
    end
    if isempty(atlastemplate) & exist(fullfile(atlas_filepath,'Brain.nii'),'file')
        atlastemplate = fullfile(atlas_filepath,'Brain.nii');
    elseif ~exist(atlastemplate,'file')
        atlastemplate = '';
    end
    if exist(atlastemplate,'file')
        atlastemplate = lab_read_mri(atlastemplate);
        ft_write_mri(fullfile(atlas_filepathout,[atlas_fileS '_brain.nii']),int16(atlastemplate.anatomy), ...
            'transform',atlastemplate.transform,'dataformat','nifti');
    end
end

if ~exist([MriFileout '_brainW.nii'],'file')
    if exist('MriFileout','var') & exist([MriFileout '_brain.hdr'],'file')
        brain = lab_read_mri([MriFileout '_brain.hdr']);
        mri2 = mri;
        mri2.anatomy(brain.anatomy==0) = 0;
        clearvars brain
        ft_write_mri([MriFileout '_brainW.nii'],int16(mri2.anatomy),'transform',mri2.transform,'dataformat','nifti');
        clearvars mri2
    else
        if exist('mri_filepath','var')
            settings.mrifile = fullfile(mri_filepath,mri_file);
            settings.SEGprobmaps = 0;
            settings.SEGcorrect = 1;
            mri = lab_segment_mri(settings);
        else
            mri = lab_segment_mri(mri);
        end
        mri.anatomy(mri.brain==0) = 0;
        mri = ft_convert_units(mri,'mm');
        ft_write_mri([MriFileout '_brainW.nii'],int16(mri.anatomy),'transform',mri.transform,'dataformat','nifti');
    end
end

[atlas_file,atlas_filepath] = lab_filename(AtlasFileout);
[mri_file,mri_filepath] = lab_filename(MriFileout);
if ~exist(fullfile(atlas_filepath,[mri_file '_' atlas_file '.hdr']),'file');
    if exist(fullfile(atlas_filepathout,[atlas_fileS '_brain.nii']),'file')
        disp('     Perform matching using normalize-function of spm8')
        eflags.smosrc = 8;
        eflags.smoref = 8;
        eflags.regtype = 'none';
        eflags.cutoff = 25;
        eflags.nits = 16;
        eflags.reg = 1;
        warning off %#ok<WNOFF>
        spm_normalise(fullfile(mri_filepath,[mri_file '_brainW.nii']), ...
            fullfile(atlas_filepathout,[atlas_fileS '_brain.nii']), ...
            [AtlasFileout '_T.mat'],'','', eflags);
        warning on %#ok<WNON>
        rflags.preserve = 0;
        if ~isfield(mri,'dim')
            mri.dim = size(mri.anatomy);
        end
        if ~isfield(mri,'originator')
            mri.originator = round(mri.dim/2);
        end
        rflags.bb =  [-mri.originator;mri.dim-mri.originator-1];
        rflags.vox = [1,1,1];
        rflags.interp = 0;
        rflags.wrap = [0,0,0];
        rflags.prefix = 'w';
        spm_write_sn([AtlasFileout '.nii'],[AtlasFileout '_T.mat'], rflags);
    end
    disp('     Perform coregister using spm8')
    estwrite.ref = {fullfile(mri_filepath,[mri_file '_brainW.nii'])};
    if exist(fullfile(atlas_filepath,['w' atlas_file '.nii']),'file')
        estwrite.source = {fullfile(atlas_filepath,['w' atlas_file '.nii'])};
    else
        estwrite.source = {fullfile(atlas_filepath,[atlas_file '.nii'])};
    end
    estwrite.other = {''};
    estwrite.eoptions.cost_fun = 'nmi';
    estwrite.eoptions.sep = [4 2];
    estwrite.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
    estwrite.eoptions.fwhm = [7 7];
    estwrite.roptions.interp = 0;
    estwrite.roptions.wrap = [0 0 0];
    estwrite.roptions.mask = 0;
    estwrite.roptions.prefix = 'r';
    out = spm_run_coreg_estwrite(estwrite);
    AtlasFileoutR = out.rfiles{end};
        
    Atlasout = ft_read_mri(AtlasFileoutR);
    Atlasout = ft_convert_units(Atlasout,'mm');
    Atlasout.anatomy = affine(Atlasout.anatomy,Atlasout.transform,[1 1 1],0,0,2);
    Atlasout.transform = mri.transform;
    Atlasout = lab_set_mricoord(Atlasout);
    lab_write_hdr(fullfile(atlas_filepath,[mri_file '_' atlas_file '.hdr']),Atlasout);
    if exist(fullfile(atlas_filepath,[atlas_file '.txt']),'file')
        copyfile(fullfile(atlas_filepath,[atlas_file '.txt']),fullfile(atlas_filepath,[mri_file '_' atlas_file '.txt']));
    end
    AtlasFileout = fullfile(atlas_filepath,[mri_file '_' atlas_file '.hdr']);
else
    disp('    read matched atlas file from previous run')
    AtlasFileout = fullfile(atlas_filepath,[mri_file '_' atlas_file '.hdr']);
end
Atlasout = lab_read_mri(AtlasFileout);

return
