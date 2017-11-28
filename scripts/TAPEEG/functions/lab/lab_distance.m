% calculates matrix of distances
%
% [distances] = lab_distance(LOCS)
%
% written by F. Hatz 2013

function distances = lab_distance(LOCS,pnt)

if isfield(LOCS,'locs') & isfield(LOCS.locs,'x') & ~isempty(LOCS.locs.x)
    numlocs = size(LOCS.locs.x,2);
    x = LOCS.locs.x;
    y = LOCS.locs.y;
    z = LOCS.locs.z;
elseif isfield(LOCS,'x') & ~isempty(LOCS.x)
    numlocs = size(LOCS.x,2);
    x = LOCS.x;
    y = LOCS.y;
    z = LOCS.z;
elseif isnumeric(LOCS) & size(LOCS,2) == 3
    numlocs = size(LOCS,1);
    x = LOCS(:,1);
    y = LOCS(:,2);
    z = LOCS(:,3);
else
    distances = [];
    return
end

if exist('pnt','var') & length(pnt) == 3
    distances = sqrt((x-pnt(1)).^2 + (y-pnt(2)).^2 + (z-pnt(3)).^2);
else
    distances = zeros(numlocs,numlocs);
    for i = 1:numlocs
        distances(i,:) = sqrt((x-x(i)).^2 + (y-y(i)).^2 + (z-z(i)).^2);
    end
end