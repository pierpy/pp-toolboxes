function dist = dist_hp_scalp(Vertices, HP, params, p, skip)
    % Initialize
    mx = params(1); my = params(2); mz = params(3); % Translation parameters
    x = params(4); y = params(5); z = params(6); % Rotation parameters

    % Rotate and translate HP points based on parameters
    Rx = [1 0 0 ; 0 cos(x) sin(x) ; 0 -sin(x) cos(x)]; % Rotation over x
    Ry = [cos(y) 0 -sin(y); 0 1 0; sin(y) 0 cos(y)]; % Rotation over y
    Rz = [cos(z) sin(z) 0; -sin(z) cos(z) 0; 0 0 1]; % Rotation over z
    HPRotated = (Rx*Ry*Rz*HP' + [mx;my;mz]*ones(1,size(HP,1)))'; % R*HP + T

    % Calculate minimum distance from each headpoint to Vertices
    dist = zeros(size(HP,1), 1); % Initialize distance vector (allocation)
    for m = 1:size(HP,1)
        % For more than 600 points, only use 1/4 of the headpoints
        if exist('skip', 'var') && skip && size(HP,1) > 600 && mod(m,4) > 0
            continue;
        end
        d = abs(bst_bsxfun(@minus, Vertices, HPRotated(m,:))); % Distances
        dist(m) = min(sum(d.^p, 2).^(1/p)); % Minimum p-norm of distances
    end
end