% Set landmarks by using template mri with landmarks
%   Template mri is normalized to the input mri and using the
%   transformation landmarks are warped to input-mri
%
% mriout = lab_match_mri(mri_file,mrilandmarks)
%
% F. Hatz 2014

function mriout = lab_match_mrilandmarks2mri(mri,mrilandmarks)

if ~exist('mrilandmarks','var')
    [tmp1,tmp2] = uigetfile('*.hdr;*.nii','Select MRI-Template with landmarks');
    mrilandmarks = fullfile(tmp2,tmp1);
    clearvars tmp1 tmp2
end
if ~exist('mri','var')
    [tmp1,tmp2] = uigetfile('*.hdr;*.nii','Select MRI-file');
    mri = fullfile(tmp2,tmp1);
    clearvars tmp1 tmp2
end

if ischar(mri)
    [mri,MRI_file] = lab_read_mri(mri);
    [~,MRI_filepath,~,MRI_fileS] = lab_filename(MRI_file);
    haspath = true;
elseif isfield(mri,'anatomy')
    MRI_filepath = pwd;
    MRI_fileS = 'MRItmp';
    haspath = false;
else
    return
end
MRI_file = fullfile(fullfile(MRI_filepath,'Coreg'),MRI_fileS);

if ischar(mrilandmarks)
    [mrilandmarks,LMRK_file] = lab_read_mri(mrilandmarks);
    [~,LMRK_filepath,~,LMRK_fileS] = lab_filename(LMRK_file);
elseif isfield(mrilandmarks,'anatomy')
    LMRK_filepath = pwd;
    LMRK_fileS = 'LMRKtmp';
else
    return
end
if haspath == true
    LMRK_file = fullfile(fullfile(MRI_filepath,'Coreg'),LMRK_fileS);
else
    LMRK_file = fullfile(LMRK_filepath,LMRK_fileS);
end

if exist('MRI_filepath','var')
    warning off %#ok<WNOFF>
    mkdir(fullfile(MRI_filepath,'Coreg'));
    warning on %#ok<WNON>
    MRI_fileout = fullfile(fullfile(MRI_filepath,'Coreg'),[MRI_fileS '_' LMRK_fileS '.hdr']);
end

if exist('MRI_fileout','var') & exist(MRI_fileout,'file') & isfield(mri,'landmarks');
    mriout = mri;
    return
end
if ~isfield(mrilandmarks,'landmarks')
    disp('    Abort: no landmarks defined');
end

ft_write_mri([MRI_file '.nii'],int16(mri.anatomy),'transform',mri.transform,'dataformat','nifti');
ft_write_mri([LMRK_file '.nii'],int16(mrilandmarks.anatomy),'transform',mrilandmarks.transform,'dataformat','nifti');

warning off %#ok<WNOFF>
disp('     Perform matching using normalize-function of spm8')
eflags.smosrc = 8;
eflags.smoref = 8;
eflags.regtype = 'none';
eflags.cutoff = 25;
eflags.nits = 16;
eflags.reg = 1;
spm_normalise([LMRK_file '.nii'],[MRI_file '.nii'],[MRI_fileout(1:end-4) '_T.mat'],'','', eflags);
warning on %#ok<WNON>

if ~isfield(mrilandmarks,'dim')
   mrilandmarks.dim = size(mrilandmarks.anatomy);
end
if ~isfield(mri,'originator')
    mrilandmarks.originator = round(mrilandmarks.dim/2);
end
rflags.preserve = 0;
rflags.bb = [-mrilandmarks.originator;mrilandmarks.dim - mrilandmarks.originator - 1];
rflags.vox = [1,1,1];
rflags.interp = 1;
rflags.wrap = [0,0,0];
rflags.prefix = 'w';
VO = spm_write_sn([MRI_file '.nii'],[MRI_fileout(1:end-4) '_T.mat'], rflags);
VOdat = affine(VO.dat,VO.mat,[1 1 1],0,0,2);
Rmri = mrilandmarks;
Rmri.anatomy = int16(VOdat);
lab_write_hdr(MRI_fileout,Rmri);

landmarks = mrilandmarks.landmarks;
for i = 1:length(landmarks);
    landmarks(i).pnt = spm_get_orig_coord(landmarks(i).pnt,[MRI_fileout(1:end-4) '_T.mat']);
end
save([fullfile(MRI_filepath,MRI_fileS) '.lmrk'],'landmarks');
mri.landmarks = landmarks;

% set output
if haspath == false
    mriout = mri;
else
    mriout = MRI_fileout;
end