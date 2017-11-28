% Helper script for lab_plot_IS
%
% [axes1,Hlight] = lab_show_brain(figure1,faces,vertices,facecolor,matrix,Valpha,Vmode,plotnr,plotslider,threshold)
%
% Written by F. Hatz Vumc Amsterdam 10/2012 / Neurology Basel 2015
% (distrubution only with permission of the author)

function [axes1,Hlight] = lab_show_brain(figure1,faces,vertices,facecolor,matrix,Valpha,Vmode,plotnr,plotslider,threshold)

if isnumeric(faces)
    faces{1} = faces;
end
if isnumeric(vertices)
    vertices{1} = vertices;
end
if isnumeric(facecolor)
    facecolor{1} = facecolor;
end
if length(faces) ~= length(vertices)
    disp('   Abort, number of faces must be equal to number of vertices')
    return
end

if ~exist('threshold','var')
    threshold = 1;
end
if ~exist('plotslider','var')
    plotslider = true;
end
if ~exist('plotnr','var')
    plotnr = 1;
end
if size(plotnr,2) == 3
    axes1 = subplot(plotnr(1,1),plotnr(1,2),plotnr(1,3));
    tmp = zeros(plotnr(1,2),plotnr(1,1));
    tmp(plotnr(1,3))=1;
    x = find(sum(tmp,2) > 0);
    y = find(sum(tmp,1) > 0);
    Vposition(1,1) = (x-1)*(1/plotnr(1,2)) + .02/plotnr(1,2);
    Vposition(1,2) = (plotnr(1,1)-y)*(1/plotnr(1,1)) + .05/plotnr(1,1);
    Vposition(1,3) = 0.95/plotnr(1,2);
    Vposition(1,4) = 0.92/plotnr(1,1);
    clearvars tmp
else
    axes1 = axes;
    Vposition = [.03 .03 .94 .94];
end
if ~exist('Vmode','var')
    Vmode = 'T';
end
if ~exist('Valpha','var')
    Valpha = ones(1,length(faces)) * 0.5;
    if length(faces) == 1
        Valpha = 1;
    elseif length(faces) == 2
        Valpha(end) = 1;
    end
    Valpha(Valpha > 1) = 1;
elseif isnumeric(Valpha) & length(Valpha) ~= length(faces)
    Valpha = ones(1,length(faces)) * Valpha(1);
    Valpha(Valpha > 1) = 1;
end
if ~exist('matrix','var')
    matrix = [];
end


set(axes1,'Visible','off','Parent',figure1,'Position',Vposition,'xtick',[],'ytick',[], ...
    'PlotBoxAspectRatio',[1 1 1],'DataAspectRatio',[1 1 1],'CameraViewAngle',6.47983385307092);

% Plot surface and volumes
for i = 1:length(faces)
    Valpha(Valpha > 1) = 1;
    patch('Parent',axes1,'Faces',faces{i},'Vertices',vertices{i},'edgecolor','none', ...
        'facealpha',Valpha(i),'FaceLighting','gouraud',...
        'FaceColor','flat','FaceVertexCData',facecolor{i});
end
material shiny

% Draw edges
if isfield(matrix,'edges') & ~isempty(matrix.edges)
    if ~isfield(matrix,'ealpha') | isempty(matrix.ealpha) | length(matrix.ealpha) ~= size(matrix.matrixedges,1)
        matrix.ealpha = ones(1,size(matrix.matrixedges,1));
    end
    matrix.ealpha(matrix.ealpha>1) = 1;
    if ~isfield(matrix,'esize')
        matrix.esize = ones(1,size(matrix.matrixedges,1));
    end
    esize = matrix.esize * 2;
    T = [];
    Tsort = [];
    for i = 1:size(matrix.matrixedges,1)
        if matrix.edges(i) > 0
            if size(esize,1) == 2 & esize(2,i) > esize(2,i)
                esize2 = esize(1,i) + (matrix.edges(i) * (esize(2,i)-esize(1,i)));
            else
                esize2 = esize(i);
            end
            loctmp(1,:) = matrix.xyz(matrix.matrixedges(i,1),:);
            loctmp(2,:) = matrix.xyz(matrix.matrixedges(i,2),:);
            alphatmp = matrix.edges(i) * matrix.ealpha(i);
            if size(matrix.ecolor,1) >= i
                ecolor = matrix.ecolor(i,:);
            else
                ecolor = matrix.ecolor(1,:);
            end
            if max(isnan(ecolor)) == 0
                T(1,end+1) = patchline(loctmp(:,1),loctmp(:,2),loctmp(:,3), ...
                    'linestyle','-','edgecolor',ecolor, ...
                    'linewidth',esize2,'edgealpha',alphatmp); %#ok<AGROW>
                Tsort(1,end+1) = matrix.edges(i); %#ok<AGROW>
            end
        end
    end
    clearvars loctmp
    if ~isempty(T)
        [~,tmp] = sort(Tsort);
        T = T(tmp);
        clearvars Tsort tmp
        if threshold < 1
            numel = length(T);
            val = round(threshold * numel);
            if val > 0
                set(T(1:val),'visible','off');
            end
        end
    end
end

% Draw nodes
if isfield(matrix,'nodes') & ~isempty(matrix.nodes)
    if ~isfield(matrix,'alpha') | isempty(matrix.alpha) | length(matrix.alpha) ~= length(matrix.nodes)
        matrix.alpha = ones(1,length(matrix.nodes));
    end
    matrix.alpha(matrix.alpha>1) = 1;
    if ~isfield(matrix,'nsize')
        matrix.nsize = ones(1,length(matrix.nodes));
    end
    sc = matrix.nodes(:) .* matrix.nsize(:) * 4;
    for i = 1:length(matrix.nodes)
        if sc(i) > 0
            if size(matrix.ncolor,1) >= i
                ncolor = matrix.ncolor(i,:);
            else
                ncolor = matrix.ncolor(1,:);
            end
            [x,y,z] = sphere(50);
            if max(isnan(ncolor)) == 0
                patch(surf2patch(x*sc(i)+matrix.xyz(i,1),y*sc(i)+matrix.xyz(i,2),z*sc(i)+matrix.xyz(i,3)), ...
                    'FaceColor',ncolor,'EdgeColor','none','FaceAlpha',matrix.alpha(i));
            end
        end
    end
    clearvars x y z sc
end

% Set camera & light position
if strcmp(Vmode,'L')
    view(axes1,[-90 0]) % left
    Hlight = light('Parent',axes1,'Position',[-1 -0.5 0]);
elseif strcmp(Vmode,'H')
    view(axes1,[0 0]) % back
    Hlight = light('Parent',axes1,'Position',[0 -1 0]);
elseif strcmp(Vmode,'R')
    view(axes1,[90 0]) % right
    Hlight = light('Parent',axes1,'Position',[1 -0.5 0]);
elseif strcmp(Vmode,'F')
    view(axes1,[180 0]) % front
    Hlight = light('Parent',axes1,'Position',[0 1 0]);
elseif strcmp(Vmode,'B')
    view(axes1,[180 -90]) % bottom
    Hlight = light('Parent',axes1,'Position',[0 0 -1]);
else
    view(axes1,[0 90]) % top
    Hlight = light('Parent',axes1,'Position',[0 0 1]);
end

if exist('T','var') & plotslider == 1
    U = get(gcf,'UserData');
    if ~isfield(U,'Hprint')
        U.Hprint = [];
    end
    U.Hprint(1,end+1) = uicontrol('Style', 'slider','Min',0,'Max',1,'Value',0,...
        'Position', [5 5 120 15],'Callback', {@adjust_threshold,T});
    set(gcf,'UserData',U);
end

return

function adjust_threshold(hObj,~,T)
   numel = length(T);
   val = round(get(hObj,'Value') * numel);
   if val > 0
       set(T(1:val),'visible','off');
   end
   if val < numel
       set(T(val+1:numel),'visible','on');
   end
return