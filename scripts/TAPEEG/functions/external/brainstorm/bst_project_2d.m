function [X,Y] = bst_project_2d(x,y,z)
% BST_PROJECT_2D: Project a set of 3D points (EEG or MEG sensors) on a 2D surface.
%
% USAGE:  [X,Y] = bst_project_2d(x,y,z);

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
% Authors: François Tadel, 2009-2010

% Spheric coordinates
z = z - max(z);
[TH,PHI,R] = cart2sph(x, y, z);
% Remove the too smal values for PHI
PHI(PHI < 0.001) = 0.001;
% Flat projection
R2 = R ./ cos(PHI) .^ .2;
[X,Y] = pol2cart(TH,R2);





