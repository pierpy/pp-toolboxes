function mri = lab_data2mri(data,mri)

if ~isfield(header,'dim') | isempty(header.dim)
    disp('   Abort: missing information about dimension')
    return
end
if size(data,1) ~= (mri.dim(1) * mri.dim(2) * mri.dim(3))
    disp('   Abort: Information about dimension not matching to data')
    return
end
if isfield(mri,'anatomy') & isempty(mri.anatomy)
    mri.anatomy = reshape(data,[mri.dim(1) mri.dim(2) mri.dim(3) size(data,2)]);
elseif isfield(mri,'img') & isempty(mri.img)
    mri.img = reshape(data,[mri.dim(1) mri.dim(2) mri.dim(3) size(data,2)]);
else
    disp('   Abort: no valid input data')
    return
end