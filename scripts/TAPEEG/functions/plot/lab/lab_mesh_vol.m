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

if ~exist('threshold','var') | ~isnumeric(threshold)  | length(threshold) ~= 1
    if length(unique(vol)) > 2
        maxval = max(vol(:));
        minval = min(vol(:));
        if round(minval) ~= minval | round(maxval) ~= maxval
            settings.threshold =  round(100*(minval+maxval) / 2) / 100;
            Prompt = {['Set Threshold (min: ' num2str(minval,'%.2f') ' max: ' num2str(maxval,'%.2f') ')'], ...
                '';'Threshold','threshold'};
        else
            settings.threshold =  (minval+maxval) / 2;
            Prompt = {['Set Threshold (min: ' num2str(minval) ' max: ' num2str(maxval) ')'], ...
                '';'Threshold','threshold'};
        end
        Formats(1,1).type = 'text';
        Formats(2,1).type = 'edit';
        Formats(2,1).format = 'float';
        Formats(2,1).limits = [-inf inf];
        settings = inputsdlg(Prompt,'Threshold',Formats,settings);
        pause(0.2);
        if isempty(settings)
            settings.threshold = (minval+maxval) / 2;
        end
        threshold = settings.threshold;
        clearvars settings Prompt Formats maxval minval
    else
        threshold = (min(vol(:))+max(vol(:))) / 2;
    end
end
vol(vol<threshold) = 0;
vol(vol>=threshold) = 1;
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
