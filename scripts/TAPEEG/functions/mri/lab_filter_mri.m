% Anisotropic filtering of mri using 'segm_aniso_filtering'
%
% [mri,settings] = lab_filter_mri(mri,settings)

function [mri,settings] = lab_filter_mri(mri,settings)

if exist('mri','var') & ischar(mri) & exist(mri,'file')
    [~,mri_filepath,~,mri_fileS] = lab_filename(mri);
    if exist(fullfile(mri_filepath,[mri_fileS '_filt.hdr']),'file')
        disp('    Use filtered mri from previous run')
        mri = fullfile(mri_filepath,[mri_fileS '_filt.hdr']);
        return
    end
    clearvars mri_filepath mri_fileS
    [~,mrifile] = lab_read_mri(mri);
    mri = load_nii(mrifile);
    img = mri.img;
elseif ~exist('mri','var') | ~isstruct(mri) | (~isfield(mri,'anatomy') & ~isfield(mri,'img'))
    [~,mrifile] = lab_read_mri;
    mri = load_nii(mrifile);
    img = mri.img;
elseif isfield(mri,'anatomy')
    img = mri.anatomy;
else
    img = mri.img;
end

if ~exist('settings','var') | isempty(settings) | ~isfield(settings,'ts')
    settings.ts = 0.0625;
    settings.iter = 5;
    settings.cond = 3;
end

disp('    Anisotropic filtering of mri')
Simg = size(img);
img = segm_aniso_filtering(img,settings.ts,settings.iter,settings.cond);
img = img(1:Simg(1),1:Simg(2),1:Simg(3));
if isfield(mri,'img')
    mri.img = img;
elseif isfield(mri,'anatomy')
    mri.anatomy = img;
end

if exist('mrifile','var')
    [~,mri_filepath,~,mrifileS] = lab_filename(mrifile);
    mrifile = fullfile(mri_filepath,[mrifileS '_Filt.hdr']);
    lab_write_hdr(mrifile,mri);
    mri = mrifile;
end