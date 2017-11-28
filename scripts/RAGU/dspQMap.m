function [c,imap,xm,ym] = dspQMap(map,ChanPos,varargin)
% dspCMap - Display topographic scalp maps
% ----------------------------------------
% Copyright 2009-2011 Thomas Koenig
% distributed under the terms of the GNU AFFERO General Public License
%
% Usage: dspCMap(map,ChanPos,CStep)
%
% map is the 1xN voltage map to display, where N is the number of electrode
%
% ChanPos contains the electrode positions, either as Nx3 xyz coordinates,
% or as structure withRadius, Theta and Phi (Brainvision convention)
%
% CStep is the size of a single step in the contourlines, a default is used
% if this is not set
%
% There are a series of options that can be set using parameter/value
% pairs:
%
% 'Resolution'    controls the resolution of the interpolation
% 'Label'         shows the electrode positions and labels them
% 'NoScale'       whether or not a scale is being shown

if nargin > 2
    if iscell(varargin{1})
        varargin = varargin{1};
    end
end

if isstruct(ChanPos)
    [x,y,z] = VAsph2cart(ChanPos);
else
    if size(ChanPos,1) == 3
        ChanPos = ChanPos';
    end
    x = -ChanPos(:,2)';
    y =  ChanPos(:,1)';
    z =  ChanPos(:,3)';
end

hold off
%cla

if (nargin < 2)
    error('Not enough input arguments');
end

CStep = [];

if (nargin > 2)
    if isnumeric(varargin{1})
        CStep = varargin{1};
    end
end    

if vararginmatch(varargin,'Resolution')
    res = varargin{vararginmatch(varargin,'Resolution')+1};
else
    res = 1;
end

if vararginmatch(varargin,'Step')
    CStep = varargin{vararginmatch(varargin,'Step')+1};
end

if vararginmatch(varargin,'Border')
    DoBorder = varargin{vararginmatch(varargin,'Border')+1};
else
    DoBorder = 0;
end


if vararginmatch(varargin,'Background')
    Background  = varargin{vararginmatch(varargin,'Background')+1};
else
    Background  = get(gcf,'Color');
end
Background = Background(:)';


r = sqrt(x.*x + y.* y);
r(r == 0) = 1;
Theta = acos(z) / pi * 180;

pxG = x./r.*Theta;
pyG = y./r.*Theta;

xmx = max(abs(pxG));
ymx = max(abs(pyG));

xa = -xmx:res:xmx;
ya = -ymx:res:ymx;

[xm,ym] = meshgrid(xa,ya);

imap = griddata(pxG,pyG,map,xm,ym,'v4');
vmap = griddata(pxG,pyG,map,xm,ym,'linear');

idx = isnan(vmap);

imap(idx) = NaN;

if DoBorder == 1
    bMat = zeros(size(vmap));
    sx = size(bMat,2);
    sy = size(bMat,1);
    rgthit = isnan(vmap(:,2:sx)) & ~isnan(vmap(:,1:sx-1));
    bMat(:,2:sx) = rgthit;
    lfthit = isnan(vmap(:,1:sx-1)) & ~isnan(vmap(:,2:sx));
    bMat(:,1:sx-1) = bMat(:,1:sx-1) + lfthit;
    uphit = isnan(vmap(2:sy,:)) & ~isnan(vmap(1:sy-1,:));
    bMat(2:sy,:) = bMat(2:sy,:) + uphit;
    dwnhit = isnan(vmap(1:sy-1,:)) & ~isnan(vmap(2:sy,:));
    bMat(1:sy-1,:) = bMat(1:sy-1,:) + dwnhit;
    bidx = bMat > 0;
end


if isempty(CStep)
    CStep = max(abs(map))/8;
end

c(1,1,:) = Background;

cmap = repmat(1,[size(imap,1) size(imap,2) 3]);
%imap = (fix(imap / CStep))/8;
imap = (ceil(abs(imap) / CStep))/8 .* sign(imap);

ipos = imap;
ipos(ipos < 0) = 0;
ineg = imap;
ineg(ineg > 0) = 0;

rmap = 1        + ineg;
gmap = 1 - ipos + ineg;
bmap = 1 - ipos;

rmap(idx) = c(1);
gmap(idx) = c(2);
bmap(idx) = c(3);


if DoBorder == 1
    rmap(bidx) = 0;
    gmap(bidx) = 0;
    bmap(bidx) = 0;
end

cmap(:,:,1) = rmap;
cmap(:,:,2) = gmap;
cmap(:,:,3) = bmap;

cmap(cmap < 0) = 0;
cmap(cmap > 1) = 1;

    
image(cmap);

hold off
axis equal
axis off
axis xy