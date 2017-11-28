% Script to read mri-files. Different file-formats will be loaded and
% converted to a nifti hdr/img format with voxel size of 1x1x1 mm, RAS
% orientation
%
% [mri,filenameout] = lab_read_mri(filename)
%
% written by F. Hatz 2012

function [mri,FilenameOut] = lab_read_mri(Filename,skiphdr)

if ~exist('skiphdr','var')
    skiphdr = 0;
end
if ~exist('Filename','var')
    [MRI_file,MRI_filepath]=uigetfile('*.hdr;*.nii;*.fif;*.fiff;*.img;*.dcm','Select MRI-file');
    Filename = fullfile(MRI_filepath,MRI_file);
end
if ~exist(Filename,'file')
    mri = [];
    FilenameOut = [];
    return
end

if skiphdr == 0
    FilenameOut = lab_convertmriname(Filename);
    [Filename,Filepath,Fileformat] = lab_filename(Filename);
else
    FilenameOut = Filename;
end
if ~exist(FilenameOut,'file')
    if strcmp(Fileformat,'dcm') | strcmp(Fileformat,'ima') | strcmp(Fileformat,'img')
        try
            mri = ft_read_mri(fullfile(Filepath,Filename),'format','dicom');
        catch %#ok<CTCH>
            mri = [];
            FilenameOut = [];
            return
        end
    else
        try
            mri = ft_read_mri(fullfile(Filepath,Filename));
        catch %#ok<CTCH>
            mri = [];
            FilenameOut = [];
            return
        end
    end
    mri = ft_convert_units(mri,'mm');
    mri = lab_set_mricoord(mri,true);
    if isa(mri.anatomy,'uint8')
        mri.anatomy = int8(mri.anatomy);
    end
    lab_write_hdr([FilenameOut(1:end-4) '.hdr'],mri);
end
mri = ft_read_mri(FilenameOut);
if isfield(mri,'hdr') & isfield(mri.hdr,'descrip') & length(mri.hdr.descrip) > 2 & ...
        strcmp(mri.hdr.descrip(end-2:end),'RAS') & mri.transform(1) ~= 1
    mri.transform(1) = 1;
    mri.transform(1,4) =  -mri.transform(1,4);
end
mri = lab_set_mricoord(mri);

[~,Filepath,~,FilenameS] = lab_filename(FilenameOut);
if exist(fullfile(Filepath,[FilenameS '.txt']),'file')
    disp(['     read labels (' FilenameS '.txt)'])
    fid = fopen(fullfile(Filepath,[FilenameS '.txt']),'r');
    if fid > 0
        labels = {};
        tline=fgetl(fid);
        while ~isnumeric(tline)
            labels{end+1,1} = tline; %#ok<AGROW>
            tline=fgetl(fid);
        end
    end
    tmp = setdiff(unique(mri.anatomy),0);
    if length(tmp) == length(labels)
        mri.labels = labels;
    end
end
if exist(fullfile(Filepath,[FilenameS '.lmrk']),'file')
    disp(['     read landmarks (' FilenameS '.lmrk)'])
    Lmrk = load(fullfile(Filepath,[FilenameS '.lmrk']),'-mat');
    if isfield(Lmrk,'landmarks')
        mri.landmarks = Lmrk.landmarks;
    end
end

return