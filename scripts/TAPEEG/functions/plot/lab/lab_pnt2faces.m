function faces = lab_pnt2faces(pnt,maxdist)

try 
    faces = projecttri(pnt,'convhull');
catch %#ok<CTCH>
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