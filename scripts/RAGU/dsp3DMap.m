function [c,imap,xm,ym] = dsp3DMap(map,ChanPos,varargin)
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

MinZ = min(z);

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

Theta = acos(z);
r = sqrt(x.*x + y.* y);
r(r == 0) = 1;

pxG = x./r.*Theta;
pyG = y./r.*Theta;

xmx = max(abs(pxG));
ymx = max(abs(pyG));

res = res / 180 * pi;

xa = -xmx:res:xmx;
ya = -ymx:res:ymx;

[xm,ym] = meshgrid(xa,ya);

imap = griddata(pxG,pyG,map,xm,ym,'v4');
vmap = griddata(pxG,pyG,map,xm,ym,'linear');
idx = isnan(vmap);

imap(idx) = NaN;

if isempty(CStep)
    CStep = max(abs(map(:)))/8;
end
if DoBorder ~=0
    nContours = floor(max(abs(map(:)))/CStep);
    (-nContours:nContours)*CStep;
    myCont = contourc(xa,ya,imap,(-nContours:nContours)*CStep);
    myCont(:,1);
else
    myCont = [];
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

cmap(:,:,1) = rmap;
cmap(:,:,2) = gmap;
cmap(:,:,3) = bmap;

cmap(cmap < 0) = 0;
cmap(cmap > 1) = 1;

rm = sqrt(xm.^2 + ym.^2);
ext = sin(rm);
zm = cos(rm);
zm(zm < MinZ) = NaN;
xm(idx) = NaN;
ym(idx) = NaN;

xm = xm * 0.8;
%zm = cos(sqrt(xm.^2 + ym.^2));

cla
hcmenu = uicontextmenu;
% Define callbacks for context menu items that change linestyle


% Define the context menu items and install their callbacks
uimenu(hcmenu, 'Label', 'Left view', 'Callback', 'view(-90,0)');
uimenu(hcmenu, 'Label', 'Right view', 'Callback', 'view(90,0)');
uimenu(hcmenu, 'Label', 'Front view', 'Callback', 'view(180,0)');
uimenu(hcmenu, 'Label', 'Back view', 'Callback', 'view(0,0)');
uimenu(hcmenu, 'Label', 'Top view', 'Callback', 'view(0,90)');

uimenu(hcmenu, 'Label', 'Left back view', 'Separator','on','Callback', 'view(-37.5,30)');
uimenu(hcmenu, 'Label', 'Right back view','Callback', 'view(37.5,30)');
uimenu(hcmenu, 'Label', 'Left front view','Callback', 'view(-127.5,30)');
uimenu(hcmenu, 'Label', 'Right front view','Callback', 'view(127.5,30)');

uimenu(hcmenu,'Label', 'Light from left', 'Separator','on','Callback', 'camlight(''Left'')');
uimenu(hcmenu,'Label', 'Light from right','Callback', 'camlight(''Right'')');
uimenu(hcmenu,'Label', 'No lights', 'Callback', 'delete(findall(gcf,''Type'',''light''))');
uimenu(hcmenu,'Label', 'Look at me', 'Separator','on','Callback',{@LookAtMe,gcf});
uimenu(hcmenu,'Label', 'Blue eyes', 'Separator','on','Callback',{@EyeColor,gcf,[0.2 0.2 1]});
uimenu(hcmenu,'Label', 'Brown eyes','Callback',{@EyeColor,gcf,[153/255 76/255 0]});


h = surface(xm./rm.*ext,ym./rm.*ext,zm,cmap,'EdgeColor','none','uicontextmenu',hcmenu);
f = 0.99;
innerc = ones(size(cmap)) * 0.2;
surface(xm./rm.*ext*f,ym./rm.*ext*f,zm*f,innerc,'EdgeColor','none');

set(h,'Tag','Scalp');

while(~isempty(myCont))
    nPts = myCont(2,1);
    cl_x = myCont(1,2:nPts);
    cl_y = myCont(2,2:nPts);
    myCont(:,1:(nPts+1)) = [];
    cl_r = sqrt(cl_x.^2 + cl_y.^2);
    cl_ext = sin(cl_r) * 1.001;
    cl_z = cos(cl_r);
    h = line(cl_x./cl_r.*cl_ext*0.8,cl_y./cl_r.*cl_ext,cl_z);
    set(h,'Color',[0 0 0],'LineWidth',DoBorder);    
end
view(3);
material dull
%camlight right

ShowEye(-0.3,0.65,-0.35,0.2,hcmenu);
ShowEye( 0.3,0.65,-0.35,0.2,hcmenu);




%item2 = uimenu(hcmenu, 'Label', 'dotted', 'Callback', hcb2);
%item3 = uimenu(hcmenu, 'Label', 'solid',  'Callback', hcb3);
% Locate line objects
% Attach the context menu to each line
set(gcf,'uicontextmenu',hcmenu);

hold off
axis equal
axis off
axis vis3d
end

function LookAtMe(hObject,eventdata,figh)

eyeh = findall(figh,'Tag','Eye');

[az,el] = view;

for i = 1:numel(eyeh)
    ud = get(eyeh(i),'Userdata');
    
    rotate(eyeh(i),[0 0 1],az-ud.dir(1),ud.pos);
    rotate(eyeh(i),[1 0 0],el-ud.dir(2),ud.pos);
    
    ud.dir(1) = az;
    ud.dir(2) = el;
    set(eyeh(i),'Userdata',ud);
end

end

function EyeColor(hObject,eventdata,figh,col)

eyeh = findall(figh,'Tag','Eye');

for i = 1:numel(eyeh)
    ud = get(eyeh(i),'Userdata');
    
    c = get(eyeh(i),'CData');
    c(ud.icol,ud.jcol,1) = col(1);
    c(ud.icol,ud.jcol,2) = col(2);
    c(ud.icol,ud.jcol,3) = col(3);
    set(eyeh(i),'CData',c);
end

end




function ShowEye(xpos,ypos,zpos,rad,contextmenu)
[x,z,y] = sphere(100);
[i1,j1] = find(y > 0.7);
[i2,j2] = find(y > 0.9);
[i3,j3] = find(y > 0.7 & y <= 0.9);
x = x*rad + xpos;
y = y*rad + ypos;
z = z*rad + zpos;
c = ones(size(x,1),size(x,2),3)-0.08;

c(i1,j1,1) = 0.2;
c(i1,j1,2) = 0.2;
c(i1,j1,3) = 1;

c(i2,j2,:) = 0;

ud.pos = [xpos ypos zpos];
ud.icol = i3;
ud.jcol = j3;
ud.dir = [180 0];
surface(x,y,z,c,'EdgeColor','none','Tag','Eye','Userdata',ud,'uicontextmenu',contextmenu);
end