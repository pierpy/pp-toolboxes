function [Transform,locs] = lab_calc_transform(p,locs)

if size(p,1) == 4 & size(p,2) == 4
    Transform = p;
elseif size(p,1) == 3 & size(p,2) == 3
    Transform = p;
    Transform(1:4,4) = zeros(4,1);
else
    mx = p(1);
    my = p(2);
    mz = p(3);
    x = p(4);
    y = p(5);
    z = p(6);
    % Rotation
    Rx = [1 0 0 ; 0 cos(x/100) sin(x/100) ; 0 -sin(x/100) cos(x/100)];
    Ry = [cos(y/100) 0 -sin(y/100); 0 1 0; sin(y/100) 0 cos(y/100)];
    Rz = [cos(z/100) sin(z/100) 0; -sin(z/100) cos(z/100) 0; 0 0 1];
    R = Rx*Ry*Rz;
    % Translation
    T = [mx; my; mz];
    
    Transform = zeros(4,4);
    Transform(1:3,1:3) = R;
    Transform(1:3,4) = T;
end

if exist('locs','var')
    if size(locs,1) == 3
        locs = Transform * [locs;ones(1,size(locs,1))];
        locs = locs(1:3,:);
    else
        locs = Transform * [locs ones(size(locs,1),1)]';
        locs = locs(1:3,:)';
    end
else
    locs = [];
end

return