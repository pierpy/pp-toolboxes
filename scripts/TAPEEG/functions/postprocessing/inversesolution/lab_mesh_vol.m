% Mesh vol, input is either mri-structure or volume (x,y,z)
%
% bnd = lab_mesh_vol(mri,nvert,threshold)
% bnd = lab_mesh_vol(mri,nvert)
% bnd = lab_mesh_vol(mri)
%
% default nvert = 3000
% default threshold = 0.5
%
% written by F. Hatz 2013

function bnd = lab_mesh_vol(mri,nvert,threshold)

disp('    mesh volume')

if ~exist('nvert','var') | ~isnumeric(nvert) | length(nvert) ~= 1
    nvert = 3000;
end

if isstruct(mri) & isfield(mri,'transform')
    transform = mri.transform;
elseif isstruct(mri) & isfield(mri,'originator')
    transform = eye(4);
    transform(1:3,4) = -mri.originator(1,1:3)';
else
    transform = eye(4);
end

if isnumeric(mri) & size(mri,3) > 1
    vol = mri;
elseif isstruct(mri) & isfield(mri,'anatomy')
    vol = mri.anatomy;
elseif isstruct(mri) & isfield(mr,'img')
    vol = mri.img;
else
    disp('Abort meshing, no valid volume information found')
    bnd = [];
    return
end

if exist('threshold','var') & isnumeric(threshold) & length(threshold) == 1
    vol(vol<threshold) = 0;
    vol(vol>=threshold) = 1;
else
    vol = (vol - min(vol(:))) / (max(vol(:))-min(vol(:)));
    vol(vol<0.5) = 0;
    vol(vol>=0.5) = 1;
end
vol = logical(vol);
vol = permute(vol,[2 1 3]);
[face,node] = isosurface(vol);
if ~isempty(node)
    [face,node] = reducepatch(face,node,10000/size(node,1));
    node = [node(:,2) -node(:,1) node(:,3)];
end
if ~isempty(node)
    node = (transform *[node ones(size(node,1),1)]')';
    bnd.pnt = node(:,1:3);
    bnd.tri = face;
else
    bnd = [];
end

return
