% Plot meshes
%
% fig1 = lab_plot_mesh(bnd,settings)
%
%       bnd                   vertices = array(number of vertices,3)
%                             volume = array(x,y,z)
%                             structure = bnd.nodes: vertices / bnd.tri: faces
%       settings.facecolor    color for faces
%       settings.plotfaces    1 = faces visible
%       settings.plotedges    1 = edges visible
%       settings.plotdots     1 = plot vertices as dots
%       settings.filename     Filename to store figure as picture
%       settings.alpha       alpha-level for mesh
%       settings.invisible    1 = figure not visible
%
% For plotting multiple meshes to single figure, call function repeatedly
%
% Written by F. Hatz 2013

function [fig1,Mhandle,Thandle,Hlight] = lab_plot_mesh(bnd,settings)

fig1 = [];
Mhandle = [];
Thandle = [];
Hlight = [];

if ~exist('settings','var') | ~isfield(settings,'facecolor')
    if ~exist('bnd','var')
        [settings,skipprocessing] = lab_set_plot_mesh;
        if skipprocessing == 1
            return
        end
        bnd = settings.bnd;
    else
        [settings,skipprocessing] = lab_set_plot_mesh(0);
        if skipprocessing == 1
            return
        end
    end
end

% correct missing settings
if isfield(settings,'plotdots') & settings.plotdots == true
    settings.plotedges = false;
    settings.plotfaces = false;
else
    settings.plotdots = false;
end
if ~isfield(settings,'dotscolor')
    settings.dotscolor = [];
end
if ~isfield(settings,'sizedots')
    settings.sizedots = 5;
end
if ~isfield(settings,'plotfaces')
    settings.plotfaces = true;
end
if ~isfield(settings,'plotedges')
    settings.plotedges = false;
end
if ~isfield(settings,'edgecolor') | isempty(settings.edgecolor) | ...
        ~isnumeric(settings.edgecolor) | length(settings.edgecolor) ~= 3
    settings.edgecolor = [0 0 0];
end
if ~isfield(settings,'alpha')
    settings.alpha = [];
end
if ~isfield(settings,'gapfactor')
    settings.gapfactor = 2;
end
if ~isfield(settings,'plotlabels')
    settings.plotlabels = false;
end
if ~isfield(settings,'labelsize')
    settings.labelsize = 20;
end
if ~isfield(settings,'labeldistance')
    settings.labeldistance = 1;
end

if ~isfield(settings,'facecolor')
    settings.facecolor = [];
elseif ischar(settings.facecolor)
    if strcmp(settings.facecolor,'y')
        settings.facecolor = [1 1 0];
    elseif strcmp(settings.facecolor,'m')
        settings.facecolor = [1 0 1];
    elseif strcmp(settings.facecolor,'c')
        settings.facecolor = [0 1 1];
    elseif strcmp(settings.facecolor,'r')
        settings.facecolor = [1 0 0];
    elseif strcmp(settings.facecolor,'g')
        settings.facecolor = [0 1 0];
    elseif strcmp(settings.facecolor,'b')
        settings.facecolor = [0 0 1];
    elseif strcmp(settings.facecolor,'w')
        settings.facecolor = [1 1 1];
    elseif strcmp(settings.facecolor,'k')
        settings.facecolor = [0 0 0];
    else
        settings.facecolor = [0 1 0];
    end
end
if ~isfield(settings,'invisible')
    settings.invisible = 0;
end
if ~isfield(settings,'filename')
    settings.filename = [];
end

if ~isstruct(bnd) & size(bnd,2) == 3
    tmp.pnt = bnd;
    if size(bnd,1) > 2000
        tmp.tri = lab_pnt2faces(bnd);
    elseif size(bnd,1) > 3
        distances = lab_distance(bnd);
        dist = sort(distances,1);
        dist = max(dist(2,:))*settings.gapfactor;
        tmp.tri = lab_pnt2faces(bnd,dist);
    else
        tmp.tri = [];
        settings.plotdots = true;
        settings.plotedges = false;
        settings.plotfaces = false;
    end
    bnd = tmp;
    clearvars tmp
elseif isnumeric(bnd) | isfield(bnd,'img') | isfield(bnd,'anatomy')
    bnd = lab_mesh_vol(bnd,10000);
end
for i = 1:size(bnd,2)
    if ~isfield(bnd(i),'pnt') | isempty(bnd(i).pnt)
        if isfield(bnd(i),'vertices')
            bnd(i).pnt = bnd(i).vertices;
        elseif isfield(bnd(i),'nodes')
            bnd(i).pnt = bnd(i).nodes;
        elseif isfield(bnd(i),'chanpos')
            bnd(i).pnt = bnd(i).chanpos;
            if isfield(bnd(i),'label')
                bnd(i).labels = bnd(i).label;
            end
        elseif isfield(bnd(i),'coilpos')
            bnd(i).pnt = bnd(i).coilpos;
            if isfield(bnd(i),'label')
                bnd(i).labels = bnd(i).label;
            end
        elseif isfield(bnd(i),'x')
            bnd(i).pnt(:,1) = bnd(i).x';
            bnd(i).pnt(:,2) = bnd(i).y';
            bnd(i).pnt(:,3) = bnd(i).z';
            bnd(i).labels = bnd(i).labels';
        else
            return
        end
    end
    if ~isfield(bnd(i),'tri') | isempty(bnd(i).tri)
        if isfield(bnd(i),'faces')
            bnd(i).tri = bnd(i).faces;
        elseif size(bnd(i).pnt,1) > 3
            if size(bnd(i).pnt,1) < 2000
                distances = lab_distance(bnd(i).pnt);
                dist = sort(distances,1);
                dist = max(dist(2,:))*settings.gapfactor;
                bnd(i).tri = lab_pnt2faces(bnd(i).pnt,dist);
            else
                bnd(i).tri = lab_pnt2faces(bnd(i).pnt);
            end
        end
    end
end

nummesh = size(bnd,2);
if isempty(settings.facecolor)
    if nummesh<= 6
        settings.facecolor = [0 1 1;1 0 1;1 1 0;1 0 0;0 1 0;0 0 1];
        settings.facecolor = settings.facecolor(1:nummesh,:);
    else
        settings.facecolor = lines(nummesh);
    end
elseif ~isnumeric(settings.facecolor)
    settings.facecolor = uisetcolor;
    settings.facecolor = repmat(settings.facecolor,nummesh,1);
elseif length(settings.facecolor) == 1
    settings.facecolor = repmat(settings.facecolor,nummesh,3);
elseif size(settings.facecolor,1) < nummesh
    tmp = nummesh - size(settings.facecolor,1);
    settings.facecolor = cat(1,settings.facecolor,repmat(settings.facecolor(end,:),tmp,1));
end
if nummesh > 1
    settings.facecolor(end,:) = [0.86 0.86 0.86];
end

if isempty(settings.dotscolor)
    if nummesh<= 6
        settings.dotscolor = [0 1 1;1 0 1;1 1 0;1 0 0;0 1 0;0 0 1];
        settings.dotscolor = settings.dotscolor(1:nummesh,:);
    else
        settings.dotscolor = lines(nummesh);
    end
elseif ~isnumeric(settings.dotscolor)
    settings.dotscolor = uisetcolor;
    settings.dotscolor = repmat(settings.dotscolor,nummesh,1);
elseif length(settings.dotscolor) == 1
    settings.dotscolor = repmat(settings.dotscolor,nummesh,3);
elseif size(settings.dotscolor,1) < nummesh
    tmp = nummesh - size(settings.dotscolor,1);
    settings.dotscolor = cat(1,settings.dotscolor,repmat(settings.dotscolor(end,:),tmp,1));
end

if isempty(settings.alpha)
    settings.alpha = ones(1,nummesh)*0.5;
elseif isnumeric(settings.alpha) & length(settings.alpha) ~= nummesh
    settings.alpha = ones(1,nummesh)*settings.alpha;
    settings.alpha(1) = 1;
end

if ~isempty(bnd) & ~isempty(bnd(1).pnt) & ~isempty(bnd(nummesh).pnt) & max(bnd(1).pnt(:,1)) >= max(bnd(nummesh).pnt(:,1))
    nummesh = nummesh:-1:1;
    settings.facecolor = settings.facecolor(nummesh,:);
    settings.alpha = settings.alpha(nummesh);
else
    nummesh = 1:nummesh;
end

plotadd = 0;
if isfield(settings,'fhandle') & ~isempty(settings.fhandle) & ishandle(settings.fhandle)
    fig1 = figure(settings.fhandle);
    axes1 = gca;
    plotadd = 1;
else
    figs = findobj('type','figure');
    for i = 1:size(figs,1)
        U = get(figs(i),'UserData');
        if isstruct(U) & isfield(U,'figuretype') & strcmp(U.figuretype,'plotmesh')
            if ~isempty(settings.filename) | settings.invisible == 1
                fig1 = figure(figs(i));
                set(fig1,'Visible','off');
            else
                fig1 = figure(figs(i));
                % set(fig1,'Color',[1 1 1]);
            end
            axes1 = gca;
            plotadd = 1;
        end
        clearvars U
    end
end

if plotadd == 0
    if ~isempty(settings.filename) | settings.invisible == 1
        fig1 = figure('Visible','off','Color',[1 1 1]);
    else
        fig1 = figure('Color',[1 1 1],'Menubar','none');
        m1 = uimenu(fig1,'Label','File');
        uimenu(m1,'Label','Save','Callback',@(~,~)lab_print_figure);
        uimenu(m1,'Label','Copy Picture','Callback','print -dmeta;');
        uimenu(m1,'Label','Close','Callback','close;');
        m2 = uimenu(fig1,'Label','Edit');
        uimenu(m2,'Label','Add Mesh','Callback',@(~,~)lab_plot_mesh);
        uimenu(m2,'Label','Black/White Background','Callback',@(~,~)set_background);
        m3 = uimenu(m2,'Label','View');
        uimenu(m3,'Label','Top','Callback',@(~,~)view_top);
        uimenu(m3,'Label','Bottom','Callback',@(~,~)view_bottom);
        uimenu(m3,'Label','Left','Callback',@(~,~)view_left);
        uimenu(m3,'Label','Right','Callback',@(~,~)view_right);
        uimenu(m3,'Label','Front','Callback',@(~,~)view_front);
        uimenu(m3,'Label','Back','Callback',@(~,~)view_back);
    end
    axes1 = axes;
    set(axes1,'Visible','off','Parent',fig1,...
        'Position',[.01 .02 .99 .95],'xtick',[],'ytick',[], ...
        'PlotBoxAspectRatio',[1 1 1],'DataAspectRatio',[1 1 1]);
    view(axes1,[-90 0]); % left
    Hlight = light('Parent',axes1,'Position',[0 0 1]);
    Hlight(2) = light('Parent',axes1,'Position',[0 0 -1]);
    U.figuretype = 'plotmesh';
    set(fig1,'Userdata',U)
    clearvars U
else
    Hlight = [];
end

Thandle = [];
Mhandle = [];
for j = nummesh
    faces = bnd(j).tri;
    vertices = bnd(j).pnt;
    
    if settings.plotfaces == true
        if settings.plotedges == true
            Mhandle(1,end+1) = patch('Parent',axes1,'Faces',faces,'Vertices',vertices,'EdgeColor',settings.edgecolor, ...
                'FaceLighting','gouraud','FaceColor','flat','facealpha',settings.alpha(j), ...
                'FaceVertexCData',repmat(settings.facecolor(j,:),size(faces,1),1)); %#ok<AGROW>
        else
            Mhandle(1,end+1) = patch('Parent',axes1,'Faces',faces,'Vertices',vertices,'EdgeColor','none', ...
                'FaceLighting','gouraud','FaceColor','flat','facealpha',settings.alpha(j), ...
                'FaceVertexCData',repmat(settings.facecolor(j,:),size(faces,1),1)); %#ok<AGROW>
        end
    end
    if settings.plotdots == true
        [x,y,z] = sphere(20);
        x = x * settings.sizedots;
        y = y * settings.sizedots;
        z = z * settings.sizedots;
        for i = 1:size(vertices,1)
            Mhandle(1,end+1) = patch(surf2patch(x+vertices(i,1),y+vertices(i,2),z+vertices(i,3)), ...
                'FaceColor',settings.dotscolor(j,:),'EdgeColor','none'); %#ok<AGROW>
        end
    end
    if settings.plotlabels == true & isfield(bnd(j),'labels')
        Tsettings.plotlabels = true;
    else
        Tsettings.plotlabels = false;
    end
    if Tsettings.plotlabels == true
        ldist = norm(mean(abs(vertices),1))*0.05;
        fsize = round(ldist*settings.labelsize)/10;
        settings.labeldistance = 1 + settings.labeldistance/10;
    end
    if Tsettings.plotlabels == true
        for i = 1:size(vertices,1)
            Thandle(i) = text(vertices(i,1)*settings.labeldistance, ...
                vertices(i,2)*settings.labeldistance, ...
                vertices(i,3)*settings.labeldistance, ...
                bnd(j).labels{i,1},'FontSize',fsize, ...
                'FontWeight','bold','FontUnit','normalized','HorizontalAlignment','center', ...
                'VerticalAlignment','middle','FontName','FixedWidth'); %#ok<AGROW>
        end
    end
    material shiny
    hold on
    rotate3d on
end
if ~isempty(settings.filename)
    view(gca,[-118 8]);
    lab_print_figure([settings.filename(1:end-4) '.jpg'],fig1);
    close;
end

   function set_background
      backgroundcolor = get(gcf,'Color');
      if backgroundcolor(1) == 1
          backgroundcolor = [0 0 0];
          textcolor = [1 1 1];
      else
          backgroundcolor = [1 1 1];
          textcolor = [0 0 0];
      end
      set(gcf,'Color',backgroundcolor);
      for T = 1:length(Thandle)
          set(Thandle(T),'Color',textcolor);
      end 
   end
   function view_top
        view(gca,[0 90]); % top
        if ~isempty(Hlight)
            set(Hlight(1),'Position',[0 0 1])
            if length(Hlight) > 1
                set(Hlight(2),'Position',[0 0 -1])
            end
        end
    end
    function view_bottom
        view(gca,[180 -90]); % bottom
        if ~isempty(Hlight)
            set(Hlight(1),'Position',[0 0 -1])
            if length(Hlight) > 1
                set(Hlight(2),'Position',[0 0 1])
            end
        end
    end
    function view_left
        view(gca,[-90 0]); % left
        if ~isempty(Hlight)
            set(Hlight(1),'Position',[-1 -0.5 0])
            if length(Hlight) > 1
                set(Hlight(2),'Position',[1 0.5 0])
            end
        end
    end
    function view_right
        view(gca,[90 0]); % right
        if ~isempty(Hlight)
            set(Hlight(1),'Position',[1 -0.5 0])
            if length(Hlight) > 1
                set(Hlight(2),'Position',[-1 0.5 0])
            end
        end
    end
    function view_front
        view(gca,[180 0]); % front
        if ~isempty(Hlight)
            set(Hlight(1),'Position',[0 1 0])
            if length(Hlight) > 1
                set(Hlight(2),'Position',[0 -1 0])
            end
        end
    end
    function view_back
        view(gca,[0 0]); % back
        if ~isempty(Hlight)
            set(Hlight(1),'Position',[0 -1 0])
            if length(Hlight) > 1
                set(Hlight(2),'Position',[0 1 0])
            end
        end
    end

end


