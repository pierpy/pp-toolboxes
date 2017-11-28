% Read mri files in nifti format
%
% mri = lab_read_nii(filename)
%
% written by F. Hatz 2012

function mri = lab_read_nii(filename)

mri = ft_read_mri(filename);

return