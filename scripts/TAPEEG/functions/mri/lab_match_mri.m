% Match mri-template to mri-file with spm8 routines
%
% [mri_file,mri] = lab_match_mri(mri_file,template_file)
%
%    mri_file        =  filename of mri-file or mri-structure (from lab_read_mri)
%    template_file   =  filename of template-file or mri-structure (from lab_read_mri)
%
% F. Hatz 2013

function [mriout_file,mriout] = lab_match_mri(mri_file,template_file,matchmode)

if ~exist('matchmode','var')
    matchmode = 'Coregister';
end

if ~exist('template_file','var') | ~exist(template_file,'file')
    [tmp1,tmp2] = uigetfile('*.hdr;*.nii;*.fif','Select Template');
    if ischar(tmp1)
        template_file = fullfile(tmp2,tmp1);
    else
        mriout_file = '';
        mriout = [];
        return
    end
    clearvars tmp1 tmp2
end
[~,template_file] = lab_read_mri(template_file);
[template_file,template_filepath,~,template_fileS] = lab_filename(template_file);

if ~exist('mri_file','var') | ~exist(mri_file,'file')
    [tmp1,tmp2] = uigetfile('*.hdr;*.nii;*.fif','Select MRI');
    if ischar(tmp1)
        mri_file = fullfile(tmp2,tmp1);
    else
        mriout_file = '';
        mriout = [];
        return
    end
    clearvars tmp1 tmp2 
end
[~,mri_file] = lab_read_mri(mri_file);
[mri_file,mri_filepath,~,mri_fileS] = lab_filename(mri_file);

mriout_filepath = fullfile(mri_filepath,'Match');
warning off %#ok<WNOFF>
mkdir(mriout_filepath);
warning on %#ok<WNON>
mriout_file = fullfile(mriout_filepath,[mri_fileS '_' template_fileS '.hdr']);

if exist(mriout_file,'file')
    mriout = lab_read_mri(mriout_file);
    return
end

if ~exist(fullfile(mriout_filepath,[mri_fileS '.nii']),'file')
    mri = lab_read_mri(fullfile(mri_filepath,mri_file));
    mri = ft_convert_units(mri,'mm');
    ft_write_mri(fullfile(mriout_filepath,[mri_fileS '.nii']),int16(mri.anatomy),'transform',mri.transform,'dataformat','nifti');
end

template = lab_read_mri(fullfile(template_filepath,template_file));
template = ft_convert_units(template,'mm');
if ~exist(fullfile(mriout_filepath,[template_fileS '.nii']),'file')
    ft_write_mri(fullfile(mriout_filepath,[template_fileS '.nii']),int16(template.anatomy),'transform',template.transform,'dataformat','nifti');
end

if strcmp(matchmode,'Coregister')
    disp('     Perform matching using coregistration-function of spm8')
    estwrite.ref = {fullfile(mriout_filepath,[template_fileS '.nii'])};
    estwrite.source = {fullfile(mriout_filepath,[mri_fileS '.nii'])};
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
    mriout_file2 = out.rfiles{1};
else
    disp('     Perform matching using normalize-function of spm8')
    eflags.smosrc = 8;
    eflags.smoref = 8;
    eflags.regtype = 'none';
    eflags.cutoff = 25;
    eflags.nits = 16;
    eflags.reg = 1;
    warning off %#ok<WNOFF>
    spm_normalise(fullfile(mriout_filepath,[template_fileS '.nii']), ...
        fullfile(mriout_filepath,[mri_fileS '.nii']), ...
        [mriout_file(1:end-4) '_T.mat'],'','', eflags);
    warning on %#ok<WNON>
    rflags.preserve = 0;
    if ~isfield(template,'dim')
        template.dim = size(template.anatomy);
    end
    if ~isfield(template,'originator')
        template.originator = round(template.dim/2);
    end
    rflags.bb =  [-template.originator;template.dim - template.originator-1];
    rflags.vox = [1,1,1];
    rflags.interp = 0;
    rflags.wrap = [0,0,0];
    rflags.prefix = 'w';
    spm_write_sn(fullfile(mriout_filepath,[mri_fileS '.nii']),[mriout_file(1:end-4) '_T.mat'],rflags);
    mriout_file2 = fullfile(mriout_filepath,['w' mri_fileS '.nii']);
end
mriout = ft_read_mri(mriout_file2);
mriout = ft_convert_units(mriout,'mm');
mriout.anatomy = affine(mriout.anatomy,mriout.transform,[1 1 1],0,0,2);
mriout.transform = mriout.transform;
mriout = lab_set_mricoord(mriout);
lab_write_hdr(mriout_file,mri);