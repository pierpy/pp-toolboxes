function [newVertices, newFaces] = tess_refine(Vertices,Faces,areaThresh,isThresh)
% TESS_REFINE: Refine a surface mesh.
% 
% USAGE:  [finalVertices, Faces] = tess_refine(Vertices, Faces, areaThresh, isThresh);
%         [finalVertices, Faces] = tess_refine(Vertices, Faces);
%
% INPUT:
%     - Vertices   : Mx3 double matrix
%     - Faces      : Nx3 double matrix
%     - areaThresh : Only split the faces that have an area above a certain threshold
%     - isThresh   : If 1, call channel_tesselate with a isThresh = 1 (remove big triangles)
%
% DESCRIPTION: Each triangle is subdivided in 4 triangles.
%     
%             /\1           
%            /  \           
%           /    \          
%          /      \        
%        4/--------\5       
%        /  \    /  \      
%       /    \6 /    \     
%     2'--------------'3   

% @=============================================================================
% This software is part of the Brainstorm software:
% http://neuroimage.usc.edu/brainstorm
% 
% Copyright (c)2000-2013 Brainstorm by the University of Southern California
% This software is distributed under the terms of the GNU General Public License
% as published by the Free Software Foundation. Further details on the GPL
% license can be found at http://www.gnu.org/copyleft/gpl.html.
% 
% FOR RESEARCH PURPOSES ONLY. THE SOFTWARE IS PROVIDED "AS IS," AND THE
% UNIVERSITY OF SOUTHERN CALIFORNIA AND ITS COLLABORATORS DO NOT MAKE ANY
% WARRANTY, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO WARRANTIES OF
% MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE, NOR DO THEY ASSUME ANY
% LIABILITY OR RESPONSIBILITY FOR THE USE OF THIS SOFTWARE.
%
% For more information type "brainstorm license" at command prompt.
% =============================================================================@
% 
% Authors: François Tadel, 2009-2011

% Parse inputs
if (nargin < 4) || isempty(isThresh)
    isThresh = 0;
end
if (nargin < 3) || isempty(areaThresh)
    areaThresh = 0;
end
% Check matrices orientation
if (size(Vertices, 2) ~= 3) || (size(Faces, 2) ~= 3)
    error('Faces and Vertices must have 3 columns (X,Y,Z).');
end

% Compute the area of all the faces
if (areaThresh > 0)
    FacesArea = zeros(1,length(Faces));
    for i = 1:length(Faces)
        A = Vertices(Faces(i,1),:);
        B = Vertices(Faces(i,2),:);
        C = Vertices(Faces(i,3),:);
        FacesArea(i) = norm(bst_cross(B-A, C-A, 2)) / 2;
    end
    % Detect the faces that have an area above the threshold
    iFacesSplit = find(FacesArea > areaThresh * std(FacesArea));
    iFacesKeep = setdiff(1:length(Faces), iFacesSplit);
else
    iFacesSplit = 1:length(Faces);
    iFacesKeep = [];
end

% ===== REFINE MESH =====
% New vertices matrix
newVertices = zeros(3*length(iFacesKeep) + 6*length(iFacesSplit), 3);
%         newFaces    = zeros(1*length(iFacesKeep) + 4*length(iFacesSplit), 3);
% Loop on faces to split
for i = 1:length(iFacesSplit)
    % Indices of available points
    P = 6*(i-1) + (1:6);
    % Copy 3 existing vertices
    Face = Faces(iFacesSplit(i),:);
    newVertices(P(1:3),:) = Vertices(Face,:);
    % Create 3 new vertices
    newVertices(P(4), :) = (Vertices(Face(1),:) + Vertices(Face(2),:)) / 2;
    newVertices(P(5), :) = (Vertices(Face(1),:) + Vertices(Face(3),:)) / 2;
    newVertices(P(6), :) = (Vertices(Face(2),:) + Vertices(Face(3),:)) / 2;
    % Create 4 new faces
%   F = 4*(i-1) + (1:4);
%   newFaces(F,:) = P([4,1,5; 4,2,6; 4,6,5; 5,6,3]);
end
% Loop on faces to keep
for i = 1:length(iFacesKeep)
    % Indices of available points
    P = 6*length(iFacesSplit) + 3*(i-1) + (1:3);
    % Copy 3 existing vertices
    newVertices(P(1:3),:) = Vertices(Faces(iFacesKeep(i),:),:);
%   % Copy one face
%   F = 4*length(iFacesSplit) + i;
%   newFaces(F,:) = P;
end

% ===== TESSELATE NEW SURFACE =====
newVertices = unique(newVertices, 'rows');
% If 3D surface, use channel_tesselate.m
if ~all(newVertices(:,3) == 0)
    newFaces = channel_tesselate(newVertices, isThresh);
% Else: flat surface, use delaunay.m
else
    newFaces = delaunay(newVertices(:,1), newVertices(:,2));
end

% ===== REMOVE DUPLICATE VERTICES =====
% [newVertices,I,J] = unique(newVertices, 'rows');
% newFaces = J(newFaces);




