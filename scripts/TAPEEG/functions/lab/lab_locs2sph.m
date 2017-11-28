function LOCS = lab_locs2sph(LOCS)

if ~exist('LOCS','var')
    LOCS = [];
end
if ~isfield(LOCS,'z')
    return
end

if min(LOCS.z) == 0 & max(LOCS.z) == 0
    [th,radius] = cart2pol(LOCS.x,LOCS.y,LOCS.z);
    LOCS.radius = radius;
    LOCS.theta = th;
else
    [th,phi,radius] = cart2sph(LOCS.x,LOCS.y,LOCS.z);
    LOCS.sph_radius = radius;
    LOCS.sph_theta = th/pi*180;
    LOCS.sph_phi = phi/pi*180;
    [~,LOCS.theta,LOCS.radius] = sph2topo([ones(length(LOCS.sph_phi),1) LOCS.sph_phi' LOCS.sph_theta'], 1, 2);
    LOCS.theta = -LOCS.theta'/180*pi;
    LOCS.radius = LOCS.radius';
end

return