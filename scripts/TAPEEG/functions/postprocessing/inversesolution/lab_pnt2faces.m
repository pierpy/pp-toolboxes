% Create faces for vertices using 'projecttri' of fieldtrip
%
% faces = lab_pnt2faces(pnt,maxdist)
%      pnt = vertices
%      maxdist = maximal distances between 2 vertices for being connected
%                by face
%
% written by F. Hatz 2012

function faces = lab_pnt2faces(pnt,maxdist)

try 
    faces = projecttri(pnt,'convhull');
catch err
    faces = [];
    return
end

if exist('maxdist','var') & ~isempty(maxdist)
   distances = zeros(size(faces));
   tmp = pnt(faces(:,1),:) - pnt(faces(:,2),:);
   distances(:,1) = sqrt(sum(tmp.^2,2));
   tmp = pnt(faces(:,1),:) - pnt(faces(:,3),:);
   distances(:,2) = sqrt(sum(tmp.^2,2));
   tmp = pnt(faces(:,2),:) - pnt(faces(:,3),:);
   distances(:,3) = sqrt(sum(tmp.^2,2));
   distances = max(distances,[],2);
   faces = faces(distances<=maxdist,:);
end