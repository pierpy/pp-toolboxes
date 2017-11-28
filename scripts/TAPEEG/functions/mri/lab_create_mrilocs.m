% Function to automatically define Nasion, Inion and Cz
%
% Not used in TAPEEG, not stable

function lab_create_mrilocs(bnd)

if ~isfield(bnd,'pnt')
    return
end
bnd = bnd(1);

% Find most anterior point in midline (should be nose)
tmp =  bnd.pnt(bnd.pnt(:,2) == max(bnd.pnt(:,2)),:);
tmp = tmp(abs(tmp(:,1)) == min(abs(tmp(:,1))),:);
nose = tmp(tmp(:,3) == min(tmp(:,3)),:);
clearvars tmp

% Find most occipital point in midline
tmp =  bnd.pnt(bnd.pnt(:,2) == min(bnd.pnt(:,2)),:);
tmp = tmp(abs(tmp(:,1)) == min(abs(tmp(:,1))),:);
occiput = tmp(tmp(:,3) == max(tmp(:,3)),:);
clearvars tmp

% Find highest point in midline
rostral = [0 0 max(bnd.pnt(:,3))];

% Create plane
plane = createPlane(nose,occiput,rostral);

% Intersect mesh and plane
disp('   intersect scalp-mesh and midline plane (please wait)')
polys = intersectPlaneMesh(plane,bnd.pnt,bnd.tri);

% Find midline contour
midline = bnd.pnt;
tmp = round(midline(:,1)*10) / 10;
midline = midline(tmp == 0,:);
clearvars tmp

% delete lowest points
midline = midline(midline(:,3) > (min(midline(:,3))+4),:);

% Sort midline
distance = lab_distance(midline);
Idx = find(midline(:,2) == max(midline(:,2)));
Idx = Idx(1);
tmp = [];
while size(tmp,1) ~= size(midline,1);
    tmp(end+1,:) = midline(Idx,:);
    distance(Idx,:) = NaN;
    Idx = find(distance(:,Idx) == min(distance(:,Idx)));
end
midline = tmp;
clearvars tmp distance Idx

% Calculate angles
for i = 2:size(midline,1)-1
    Mangle(i,1) = abs(anglePoints3d(midline(i-1:i+1,:)));
    Mangle(i,2) = midline(i,2) / (0.5 * (midline(i-1,2) + midline(i-1,2)));
end