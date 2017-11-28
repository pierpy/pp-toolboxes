% Find individual file for calculation of inverse solutions%
%   used by lab_inversesolution / lab_inversesolution_fft
%
% written by F. Hatz 2012

function [IS_file,SPI_file,ROIS_file,MRI_file] = lab_find_individualIS(IS_file,SPI_file,ROIS_file,MRI_file,Path)

if strcmp(IS_file,'nii') | strcmp(IS_file,'hdr') | strcmp(IS_file,'fif') | strcmp(IS_file,'fiff')
    domri = 1;
else
    domri = 0;
end

if strcmp(IS_file,'dicom')
    findis = dir(fullfile(fullfile(Path,'DICOM'),'*'));
    if size(findis,1) > 2
        findis = findis(3,1);
    else
        findis = dir(fullfile(fullfile(Path,'dicom'),'*'));
        if size(findis,1) > 2
            findis = findis(3,1);
        else
            findis = [];
        end
    end
else
    findis = dir(fullfile(Path,['*' IS_file]));
end
if size(findis,1) > 0
    IS_file = fullfile(Path,findis(1,1).name);
else
    IS_file = '';
end

if strcmp(SPI_file,'dicom')
    findspi = dir(fullfile(fullfile(Path,'DICOM'),'*'));
    if size(findspi,1) > 2
        findspi = findis(3,1);
    else
        findspi = dir(fullfile(fullfile(Path,'dicom'),'*'));
        if size(findspi,1) > 2
            findspi = findis(3,1);
        else
            findspi = [];
        end
    end
elseif ~isempty(SPI_file) & ~strcmp(SPI_file,' ')
    findspi = dir(fullfile(Path,['*' SPI_file]));
else
    findspi = [];
end
if size(findspi,1) > 0
    SPI_file = fullfile(Path,findspi(1,1).name);
else
    SPI_file = '';
end

if length(ROIS_file) < 5
    if ~isempty(ROIS_file) & ~strcmp(ROIS_file,' ')
        findrois = dir(fullfile(Path,['*' ROIS_file]));
    else
        findrois = [];
    end
    if size(findrois,1) > 0
        ROIS_file = fullfile(Path,findrois(1,1).name);
    else
        ROIS_file = '';
    end
    if strcmp(IS_file,ROIS_file)
        ROIS_file = '';
    end
end

if domri == 1
    [~,MRI_path,~,MRI_file] = lab_filename(IS_file);
    MRI_file = fullfile(MRI_path,[MRI_file '.hdr']);
else
    findmri = dir(fullfile(Path,['*' MRI_file]));
    if size(findmri,1) > 0
        MRI_file = fullfile(Path,findmri(1,1).name);
    else
        MRI_file = '';
    end
end

return