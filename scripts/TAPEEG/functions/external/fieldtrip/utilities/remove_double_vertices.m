function [pntR, triR] = remove_double_vertices(pnt, tri)

% REMOVE_VERTICES removes specified vertices from a triangular mesh
% renumbering the vertex-indices for the triangles and removing all
% triangles with one of the specified vertices.
%
% Use as
%   [pnt, tri] = remove_double_vertices(pnt, tri)

pnt1 = unique(pnt, 'rows');
keeppnt   = find(ismember(pnt1,pnt,'rows'));
removepnt = setdiff([1:size(pnt,1)],keeppnt);

npnt = size(pnt,1);
ntri = size(tri,1);

if all(removepnt==0 | removepnt==1)
  removepnt = find(removepnt);
end

% remove the vertices and determine the new numbering (indices) in numb
keeppnt = setdiff(1:npnt, removepnt);
numb    = zeros(1,npnt);
numb(keeppnt) = 1:length(keeppnt);

% look for triangles referring to removed vertices
removetri = false(ntri,1);
removetri(ismember(tri(:,1), removepnt)) = true;
removetri(ismember(tri(:,2), removepnt)) = true;
removetri(ismember(tri(:,3), removepnt)) = true;

% remove the vertices and triangles
pntR = pnt(keeppnt, :);
triR = tri(~removetri,:);

% renumber the vertex indices for the triangles
triR = numb(triR);