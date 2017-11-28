% segm_aniso_filtering() - Performs anisotropic filtering.
%
% Usage:
%   >> [b] = segm_aniso_filtering(out, iter, ts, cond);
%
% Inputs:
%   out - input image
%   iter - number of iterations
%   ts - sampling rate
%   cond - diffusion parameter
%
% Outputs:
%   b - output image
%
% Author: Zeynep Akalin Acar, SCCN, 2008

% Copyright (C) 2007 Zeynep Akalin Acar, SCCN, zeynep@sccn.ucsd.edu
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

function [mri] = segm_aniso_filtering(mri,ts,iter,cond)

% add slices of zero to the bottom and sides of the head
[K,L,M] = size(mri);
mri2 = zeros(K,L,M);
mri2(:,2:L+1,:) = mri;
ek = zeros(K,M);
mri2(:,1,:) = ek;
mri = mri2;
clear mri2
[K,L,M] = size(mri);
mri2 = zeros(K,L,M+2);
mri2(:,:,2:M+1) = mri;
mri = mri2;
clear mri2

% filtering
mri = matitk('FCA',[iter ts cond],mri);

