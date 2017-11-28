% Plot electrodes and mappings
%
% [settings,skipprocessing] = lab_plot_chans(data,settings,plotsingle,plotslider,plothead)
%
% data = array with values (Mappings x 1 columns)
% settings.Mappings   = filepath to Mappings.xls or Mappings-structure
% settings.LOCS       = filepath to LOCS or LOCS-structure
% settings.PLOT_file  = Filepath to store result as image
%
% PLOT.MinValue       = minimal value to plot
% PLOT.MaxValue       = maximal value to plot
% PLOT.Color          = Color for nodes
% PLOT.Size           = Size of Nodes
% PLOT.ColorE         = Color of connectivities
% PLOT.SizeE          = Size of connectivities
% PLOT.ThresholdE     = Threshold for connectivity plotting
% PLOT.Tlegend        = Thresholds for legend plotting
% PLOT.TlegendE       = Thresholds for legend plotting (connections)
% PLOT.Legend         = Legend plotting
% PLOT.LegendE        = Legend plotting (connections)
%
% plotsingle          = if = 1, allways single dots are plotted
% plotslider          = if = 1, for matrices a slider to control amount of
%                               connections is plotted
% plothead            = if = 0, no headshape is plotted
%                       if = 1, plot headshape
%                       if = 2, plot headshape and restrict surface-plotting to header
%
% code for plotting headshape, nose and ears: eeglab
%
% Written by F. Hatz 2014

function [settings,skipprocessing] = lab_plot_chans(data,PLOT,settings,plotsingle,plotslider,plothead,Legend)

if ~exist('settings','var')
    settings = [];
end
if ~exist('plotsingle','var') | isempty(plotsingle)
    plotsingle = 0;
end
if ~exist('plotslider','var') | isempty(plotsingle)
    plotslider = 0;
end
if ~exist('plothead','var') | isempty(plothead)
    plothead = 1;
end
if ~exist('Legend','var')
    Legend = {};
end

SetTmp.LOCS = settings.LOCS;
if size(data,1) == size(data,2)
    SetTmp.DATA.data = lab_extract_tril_wodiag(data);
else
    SetTmp.DATA.data = data;
end
SetTmp.DATA.data = SetTmp.DATA.data(:)';
SetTmp.PLOT = PLOT;
if isfield(settings,'Mappings')
    SetTmp.Mappings = settings.Mappings;
end
[SetTmp,skipprocessing] = lab_plot_match_numchans(SetTmp);
if skipprocessing == 1
    return
else
    settings.LOCS = SetTmp.LOCS;
    if isfield(SetTmp,'Mappings')
        settings.Mappings = SetTmp.Mappings;
    end
end
clearvars SetTmp

if size(data,1) == size(data,2)
    matrix = data;
    data = zeros(size(matrix,1),size(matrix,3));
    for i = 1:size(matrix,3)
        data(:,i) = diag(matrix(:,:,i));
    end
    plotsingle = 1;
    plotslider = 1;
else
    matrix = [];
end
if size(data,2) > 1 & size(data,2) > size(data,1)
    data = data';
end
if size(data,2) > 1
    data = mean(data,2);
end
if size(matrix,3) > 1
    matrix = mean(matrix,3);
end
Nchans = size(settings.LOCS.x,2);
if isempty(data)
    data = zeros(Nchans,1);
    settings.nodiag = 1;
end
if isfield(PLOT,'Size') & ~isempty(PLOT.Size)
    if length(PLOT.Size) == Nchans
        Nsize = PLOT.Size;
    else
        Nsize = PLOT.Size(1) * ones(1,Nchans);
    end
else
    Nsize = ones(1,Nchans);
end
if isfield(PLOT,'SizeE') & ~isempty(PLOT.SizeE) & ~isempty(matrix)
    if length(PLOT.SizeE) == (Nchans^2/2 - Nchans/2)
        SizeE = PLOT.SizeE;
    else
        SizeE = PLOT.SizeE(1) * ones(1,Nchans^2/2 - Nchans/2);
    end
else
    SizeE = Nsize(1) * ones(1,Nchans^2/2 - Nchans/2);
end
if isfield(PLOT,'Color')
    Color = PLOT.Color;
    if ~isnumeric(Color) | size(Color,1) == 1
        Color = lab_create_cmap(Color);
    end
else
    Color = lab_create_cmap([1 0 0]);
end
if isfield(PLOT,'ColorE')
    ColorE = PLOT.ColorE;
    if ~isnumeric(ColorE) | size(ColorE,1) == 1
        ColorE = lab_create_cmap(ColorE);
    end
else
    ColorE = Color;
end
if isfield(PLOT,'MinValue')
    MinValue = PLOT.MinValue;
    MaxValue = PLOT.MaxValue;
else
    MinValue = min(data(:));
    MaxValue = max(data(:));
end
if MinValue == MaxValue
    MinValue = 0;
end
if isfield(PLOT,'MinValueE')
    MinValueE = PLOT.MinValueE;
    MaxValueE = PLOT.MaxValueE;
else
    MinValueE = MinValue;
    MaxValueE = MaxValue;
end
if MinValueE == MaxValueE
    MinValueE = 0;
end

if Nchans == size(data,1)
    for i = 1:size(data,1)
        settings.Mappings.mappings{1,i} = i;
        settings.Mappings.mappingstitle{i,1} = settings.LOCS.labels{1,i};
    end
    settings.Mappings.mappingsChannels = size(data,1);
end
if ~isfield(settings,'Mappings')
    disp('Abort: Error')
    return
end
if ~isfield(settings.LOCS,'radius')
    disp('Abort: Error LOCS')
    return
end
if Nchans < settings.Mappings.mappingsChannels
    disp('Abort: Mismatch LOCS and mappings')
    return
end

LOCS = settings.LOCS;
Mappings = settings.Mappings;
Nsize = Nsize * 0.016;

% create figure if necessary
if ~isfield(PLOT,'AddPlot') | PLOT.AddPlot == 0
    U.figuretype = 'plot2D';
    U.numelectrodes = size(LOCS.x,2);
    if ~isfield(settings,'handleF') | isempty(settings.handleF) | ~ishandle(settings.handleF)
        settings.handleF = figure('Color',[1 1 1],'Renderer','zbuffer','UserData',U,'Menubar','none');
    else
        if ~isfield(get(settings.handleF,'UserData'),'figuretype');
            set(settings.handleF,'UserData',U);
        end
    end
    clearvars U
    if ~isfield(settings,'NoMenu') | settings.NoMenu == false
        m1 = uimenu(gcf,'Label','File');
        uimenu(m1,'Label','Save Picture','Callback',@(~,~)lab_print_figure);
        uimenu(m1,'Label','Copy Picture','Callback','print -dmeta;');
        uimenu(m1,'Label','Close','Callback','close;');
        if isfield(PLOT,'Legend')
            settings.handleM = uimenu(gcf,'Label','Edit');
            uimenu(settings.handleM,'Label','Toggle Legend','Callback',@(~,~)toggle_legend);
            uimenu(settings.handleM,'Label','Legend-Table','Callback',@(~,~)show_legend);
            uimenu(settings.handleM,'Label','Edit Legend','Callback',@(~,~)toggle_legend2);
        end
    end
elseif isfield(settings,'handleF') & ishandle(settings.handleF)
    set(0,'CurrentFigure',settings.handleF)
    hold on
else
    settings.handleF = gcf;
    U.figuretype = 'plot2D';
    U.numelectrodes = size(LOCS.x,2);
    set(settings.handleF,'UserData',U);
    hold on
end

% set colormap
NumCmap = size(Color,1);

% get radius for plotting
if ~isfield(LOCS,'radius')
    LOCS = lab_locs2sph(LOCS);
end
if max(abs(LOCS.radius)) > 1
    LOCS.radius = LOCS.radius / max(abs(LOCS.radius));
end
[x,y] = pol2cart(LOCS.theta,LOCS.radius);
rd = LOCS.radius;

% prepare plot
circgrid = 201;
headcolor = [0.3 0.3 0.3];
hlinewidth = 1.7;
rmax = 0.5;
hwidth = .004;
plotrad = min(1.0,max(rd)*1.05);
plotrad = max(plotrad,0.5);
intrad = min(1.0,max(rd)*1.05);
intchans = find(rd <= intrad);
if plotrad >= rmax
    headrad = rmax;  % anatomically correct
else % if plotrad < rmax
    headrad = 0;    % don't plot head
    fprintf('topoplot(): not plotting cartoon head since plotrad (%5.4g) < 0.5\n',plotrad);
end
intx  = x(intchans);
inty  = y(intchans);
squeezefac = rmax/plotrad;
% to plot all inside the head cartoon
intx = intx*squeezefac;
inty = inty*squeezefac;
x    = x*squeezefac;
y    = y*squeezefac;
hin  = squeezefac*headrad*(1- hwidth/2);
rotate = 3*pi/2;
allcoords = (inty + intx*sqrt(-1))*exp(sqrt(-1)*rotate);
intx = imag(allcoords);
inty = real(allcoords);
allcoords = (y + x*sqrt(-1))*exp(sqrt(-1)*rotate);
x = imag(allcoords);
y = real(allcoords);
circ = linspace(0,2*pi,circgrid);
rx = sin(circ);
ry = cos(circ);
headx = [[rx(:)' rx(1) ]*(hin+hwidth)  [rx(:)' rx(1)]*hin];
heady = [[ry(:)' ry(1) ]*(hin+hwidth)  [ry(:)' ry(1)]*hin];
headx2 = [rx(:)' rx(1) ]*(hin+hwidth); %#ok<NASGU>
heady2 = [ry(:)' ry(1) ]*(hin+hwidth); %#ok<NASGU>
base  = rmax-.0046;
basex = 0.18*rmax; % nose width
tip   = 1.15*rmax;
tiphw = .04*rmax; % nose tip half width
tipr  = .01*rmax; % nose tip rounding
q = .04; % ear lengthening
earx  = [.497-.005  .510  .518  .5299 .5419  .54    .547   .532   .510   .489-.005]; % rmax = 0.5
eary  = [q+.0555 q+.0775 q+.0783 q+.0746 q+.0555 -.0055 -.0932 -.1313 -.1384 -.1199];
sf    = headrad/plotrad;
X = y;
Y = -x;

% set axes
if ~isfield(PLOT,'AddPlot') | PLOT.AddPlot == 0
    set(gca,'Xlim',[-rmax rmax]*1,'Ylim',[-rmax rmax]*1,'FontSize',9,'FontUnit', ...
        'normalized','visible','off','Position',[0 0 1 1]);
    view(gca,[0 90]) % top
end

% plot head, nose and ears
mapvars = size(Mappings.mappings,2);
if mapvars == size(LOCS.x,2)
    if plotsingle == 0
        dephtvar = 0.001;
    else
        dephtvar = -0.0015;
    end
else
    dephtvar = 0.001;
end
if (~isfield(PLOT,'AddPlot') | PLOT.AddPlot == 0) & plothead > 0
    patch(headx,heady,dephtvar*ones(size(headx)),headcolor,'edgecolor',headcolor);
    % plot3(headx2,heady2,dephtvar*ones(size(headx2)),'color',headcolor,'linewidth',hlinewidth);
    hold on
    plot3([basex;tiphw;0;-tiphw;-basex]*sf,[base;tip-tipr;tip;tip-tipr;base]*sf,...
        dephtvar*ones(size([basex;tiphw;0;-tiphw;-basex])),...
        'color',headcolor,'linewidth',hlinewidth); % plot nose
    plot3(earx*sf,eary*sf,dephtvar*ones(size(earx)),'color',headcolor,'linewidth',hlinewidth); % plot left ear
    plot3(-earx*sf,eary*sf,dephtvar*ones(size(eary)),'color',headcolor,'linewidth',hlinewidth); % plot right ear
end

% plot data
if ~isfield(settings,'handleL')
    settings.handleL = NaN(Nchans,1);
end
if ~isfield(settings,'handleT')
    settings.handleT = NaN(Nchans,1);
end
savedots = [];
for i = 1:mapvars
    sel = Mappings.mappings{1,i};
    if length(sel) == 1
        savedots = [savedots i]; %#ok<AGROW>
    else
        mx = X(sel);
        my = Y(sel);
        mxc = mean(mx);
        myc = mean(my);
        mx = mx-mxc;
        my = my-myc;
        cnv = convhull(mx,my);
        mxh = mx(cnv);
        myh = my(cnv);
        mxh = mxh + mxc;
        myh = myh + myc;
        [geom,~,cpmo] = polygeom(mxh,myh);
        mxc = geom(2);
        myc = geom(3);
        if cpmo(1) > cpmo(3)
            tangle = cpmo(4) / (pi/180);
        else
            tangle = cpmo(2) / (pi/180);
        end
        clearvars geom cpmo
        mxt = mxh*10000;
        myt = myh*10000;
        mxt = round(mxt);
        myt = round(myt);
        if length(mxt) == 2 | length(unique(mxt)) == 1 | length(myt) == 2 | length(unique(myt)) == 1
            mxh = [max(mxh)+intrad/100 min(mxh)-intrad/100 min(mxh)-intrad/100 max(mxh)+intrad/100 max(mxh)+intrad/100];
            myh = [min(myh)-intrad/100 min(myh)-intrad/100 max(myh)+intrad/100 max(myh)+intrad/100 min(myh)-intrad/100];
        end
        clearvars mx my mxt myt cnv
        if ~isempty(matrix)
            xmatrix(i) = mxc; %#ok<AGROW>
            ymatrix(i) = myc; %#ok<AGROW>
        end
        if (min(data) == 0 & max(data) == 0) | (isfield(Mappings,'plotdots') & Mappings.plotdots == true)
            th = 0:pi/180:2*pi;
            mxh = mean(Nsize(sel)) * cos(th) + mxc;
            myh = mean(Nsize(sel)) * sin(th) + myc;
            tangle = 0;
            myc = myc - mean(Nsize(sel))*2;
        end
        if ~isnan(data(i,1))
            intensity = ceil((data(i,1) - MinValue)/(MaxValue - MinValue) * NumCmap);
            if intensity > NumCmap
                intensity = Color(NumCmap,:);
            elseif intensity < 1 | isnan(intensity)
                intensity = [1 1 1];
            else
                intensity = Color(intensity,:);
            end
            if strcmp(intensity,'w') | min(intensity == [1 1 1])
                if ~ishandle(settings.handleL(sel(1)))
                    settings.handleL(sel) = plot(mxh,myh,'k-');
                end
            else
                settings.handleL(sel) = fill(mxh,myh,intensity,'EdgeColor','k');
            end
        end
        if ~isempty(matrix)
            tangle = 0;
        elseif tangle < -90
            tangle = tangle + 180;
        elseif tangle > 90
            tangle = tangle - 180;
        end
        if ~ishandle(settings.handleT(sel(1)))
            fsize = round(mean(Nsize(sel))*2500)/8;
            if isfield(Mappings,'mappingstitleS') & isfield(Mappings,'shortnames') & Mappings.shortnames == true
                fsize = round(mean(Nsize(sel))*2500)/3;
                settings.handleT(sel) = text(mxc,myc,0.002,Mappings.mappingstitleS{i,1},'FontSize',fsize, ...
                    'FontUnit','normalized','HorizontalAlignment','center', ...
                    'VerticalAlignment','middle','Rotation',tangle);
            else
                settings.handleT(sel) = text(mxc,myc,0.002,Mappings.mappingstitle{i,1},'FontSize',fsize, ...
                    'FontUnit','normalized','HorizontalAlignment','center', ...
                    'VerticalAlignment','middle','Rotation',tangle);
            end
            % g1 = get(t,'Extent');
            % set(t,'HorizontalAlignment','left','VerticalAlignment','bottom');
            % g2 = get(t,'Extent');
            % set(t,'Position',[mxc-g2(1)+g1(1) myc-g2(2)+g1(2) 0.002]);
        else
            set(settings.handleT(sel),'Position',[mxc myc 0.002], ...
                'HorizontalAlignment','center','VerticalAlignment','middle','Rotation',tangle);
        end
        clearvars mxh myh sel mxc myc
    end
end

if length(savedots) == mapvars & plotsingle == 0
    xmin = -intrad*squeezefac;
    xmax = intrad*squeezefac;
    ymin = -intrad*squeezefac;
    ymax = intrad*squeezefac;
    GRID_SCALE = 1000;
    xi = linspace(xmin,xmax,GRID_SCALE);   % x-axis description (row vector)
    yi = linspace(ymin,ymax,GRID_SCALE);
    intdata = data(intchans,1)';
    unsh = (GRID_SCALE+1)/GRID_SCALE;
    mx = intx;
    my = inty;
    mxc = mean(mx);
    myc = mean(my);
    mx = mx-mxc;
    my = my-myc;
    cnv=convhull(mx,my);
    mx = intx(cnv)*1.03;
    my = inty(cnv)*1.03;
    mdata = intdata(cnv);
    warning off %#ok<WNOFF>
    F = TriScatteredInterp([inty my]',[intx mx]',[intdata mdata]');
    warning on %#ok<WNON>
    [Xi,Yi] = meshgrid(xi,yi);
    Xi = real(Xi);
    Yi = real(Yi);
    Zi = real(F(Xi,-Yi));
    Zi = ceil((Zi - min(Zi(:)))/(max(Zi(:)) - min(Zi(:))) * NumCmap);
    Zi(Zi>NumCmap) = NumCmap;
    Zi(Zi<1) = 1;
    Idx = find(~isnan(Zi));
    Cidx = Zi(Idx);
    for i = 1:3
        tmp = NaN(size(Xi));
        tmp(Idx) = Color(Cidx,i);
        if plothead == 2
            tmp(sqrt(Xi.^2 + Yi.^2) > max(headx)) = NaN;
        end
        Zi(:,:,i) = tmp;
    end
    Xi = Xi*unsh;
    Yi = Yi*unsh;
    
    l = surface(Xi,Yi,zeros(size(Yi)),Zi,'EdgeColor','none','FaceColor','flat');
    settings.handleL = ones(Nchans,1)*l;
    settings.SURF.Xi = Xi;
    settings.SURF.Yi = Yi;
    settings.SURF.Zi = Zi;
    skipdots = 1;
    clearvars xmin xmax ymin ymax xi yi intdata Xi Yi Zi unsh mx my mxc myc mxh myh cnv tmp
else
    skipdots = 0;
end

% plot single dots
for i = savedots
    sel = Mappings.mappings{1,i};
    th = 0:pi/180:2*pi;
    mxh = mean(Nsize(sel)) * cos(th) + X(sel);
    myh = mean(Nsize(sel)) * sin(th) + Y(sel);
    if ~isempty(matrix)
        xmatrix(i) = X(sel); %#ok<AGROW>
        ymatrix(i) = Y(sel); %#ok<AGROW>
    end
    if ~isnan(data(i,1))
        intensity = ceil((data(i,1) - MinValue)/(MaxValue - MinValue) * NumCmap);
        if intensity > NumCmap
            intensity = Color(NumCmap,:);
        elseif intensity < 1 | isnan(intensity)
            intensity = [1 1 1];
        else
            intensity = Color(intensity,:);
        end
        if skipdots == 0
            if min(intensity == [1 1 1])
                if ~ishandle(settings.handleL(i))
                    settings.handleL(i) = fill(mxh,myh,intensity,'EdgeColor','k');
                end
            elseif ~isempty(matrix)
                settings.handleL(i) = fill(mxh,myh,intensity,'EdgeColor','k');
            else
                settings.handleL(i) = fill(mxh,myh,intensity,'EdgeColor',intensity);
            end
        end
    end
    if ~ishandle(settings.handleT(i))
        if skipdots == 0
            fsize = round(mean(Nsize(sel))*3350)/10;
        else
            fsize = round(mean(Nsize(sel))*3700)/10;
        end
        if fsize > 0
            t = text(X(i),Y(i),0.002,Mappings.mappingstitle{i,1},'FontSize',fsize, ...
                'FontWeight','bold','FontUnit','normalized','HorizontalAlignment','center', ...
                'VerticalAlignment','middle','FontName','FixedWidth');
            % g1 = get(t,'Extent');
            % set(t,'HorizontalAlignment','left','VerticalAlignment','bottom');
            % g2 = get(t,'Extent');
            % set(t,'Position',[X(i)-g2(1)+g1(1) Y(i)-g2(2)+g1(2) 0.002]);
            settings.handleT(i) = t;
        end
    end
    clearvars myh mxh th
end

% plot connections
if ~isempty(matrix) & size(matrix,1) == size(data,1)
    MatrixMin = min(matrix(:));
    MatrixMax = max(matrix(:));
    [matrixline,~,matrix_x,matrix_y] = lab_extract_tril_wodiag(matrix);
    [matrixline,tmp] = sort(matrixline);
    matrix_x = matrix_x(tmp);
    matrix_y = matrix_y(tmp);
    matrixline = [matrixline matrix_x matrix_y];
    tmp2 = matrixline(:,1) ~= 0;
    SizeE = (1.55 - 0.002*size(matrix,1)) * SizeE;
    if ~isempty(tmp2)
        matrixline = matrixline(tmp2,:);
        NumCon = size(matrixline,1);
        settings.handleE = zeros(1,NumCon);
        for i=1:NumCon
            intensity = ceil((matrixline(i,1) - MinValueE)/(MaxValueE - MinValueE) * NumCmap);
            if intensity > NumCmap
                intensity = ColorE(NumCmap,:);
            elseif intensity < 1 | isnan(intensity)
                intensity = [1 1 1];
            else
                intensity = ColorE(intensity,:);
            end
            Z = [dephtvar dephtvar] + 0.001 * (i/(NumCon+1));
            if min(intensity) < 1
                settings.handleE(1,i) = line([xmatrix(matrixline(i,2)) xmatrix(matrixline(i,3))], ...
                    [ymatrix(matrixline(i,2)) ymatrix(matrixline(i,3))],Z, ...
                    'LineWidth',SizeE(i),'Color',intensity);
            else
                settings.handleE(1,i) = line([xmatrix(matrixline(i,2)) xmatrix(matrixline(i,3))], ...
                    [ymatrix(matrixline(i,2)) ymatrix(matrixline(i,3))],Z, ...
                    'LineWidth',SizeE(i),'Color',[1 1 1],'Visible','off');
            end
        end
        if isfield(settings,'handleM') & ishandle(settings.handleM)
            uimenu(settings.handleM,'Label','Set Threshold','Callback',@(~,~)adjust_threshold_menu);
        end
    else
        settings.handleE = [];
    end
else
    settings.handleE = [];
end

% plot legend
U = get(gcf,'UserData');
if ~isempty(Legend)
    U.Legend = Legend;
elseif isfield(PLOT,'Legend') & ~isempty(PLOT.Legend)
    if ~isfield(U,'Legend')
        U.Legend{1} = PLOT.Legend;
    else
        U.Legend{end+1} = PLOT.Legend;
    end
elseif ~isfield(U,'Legend')
    U.Legend = {};
end
if ~isempty(U.Legend)
    if isfield(U,'Lhandle') & ~isempty(U.Lhandle) & ishandle(U.Lhandle)
        Lflag = true;
    else
        Lflag = false;
    end
    for i = 1:length(U.Legend)
        for i2 = 1:length(U.Legend{i})
            if isfield(U.Legend{i}(i2),'Mode') & ~strcmp(U.Legend{i}(i2).Mode,'Colorbar')
                Lflag = true;
            end
        end
    end
    if plotsingle == 0 | (~isempty(matrix) & min(data) ~= max(data))
        Lflag = true;
    end
    if Lflag == true
        for M = 1:length(U.Legend)
            U = lab_plot_legend(U.Legend{M},U);
        end
    end
end
set(gcf,'UserData',U);

% set threshold of connectivities
if isfield(PLOT,'ThresholdE') & ~isempty(PLOT.ThresholdE)
    numel = length(settings.handleE);
    val = round((PLOT.ThresholdE/100) * numel);
    if val > 0
        set(settings.handleE(1:val),'visible','off');
    end
    if val < numel
        set(settings.handleE(val+1:numel),'visible','on');
    end
end

% save result to image file
if isfield(PLOT,'Title')
    set(gcf,'Name',PLOT.Title,'NumberTitle','off');
end
if isfield(settings,'PLOT_file') & ~isempty(settings.PLOT_file) & ~isnumeric(settings.PLOT_file)
    lab_print_figure(settings.PLOT_file,gcf);
end
if isfield(settings,'close') & settings.close == 1
    close;
end
if isfield(settings,'handleE') & exist('plotslider','var') & plotslider == 1
    U = get(gcf,'UserData');
    if ~isfield(U,'Hprint')
        U.Hprint = [];
    end
    Hslider = uicontrol('Style', 'slider','Min',0,'Max',1,'Value',0,...
        'Position', [5 5 120 15],'Callback', @(~,~)adjust_threshold);
    U.Hprint(1,end+1) = Hslider;
    set(gcf,'UserData',U);
    clearvars U
    if isfield(PLOT,'ThresholdE') & ~isempty(PLOT.ThresholdE)
        set(Hslider,'Value',PLOT.ThresholdE/100);
    end
else
    Hslider = [];
end

     function toggle_legend
        U2 = get(gcf,'UserData');
        if Lflag == true
            U2 = lab_plot_legend('off',U2);
            Lflag = false;
        else
            for m = 1:length(U2.Legend)
                U2 = lab_plot_legend(U2.Legend{m},U2);
            end
            Lflag = true;
        end
        set(gcf,'UserData',U2);
     end
    
     function toggle_legend2
        U2 = get(gcf,'UserData');
        if Lflag == true
            U2 = lab_plot_legend('off',U2);
            Lflag = false;
        else
            for m = 1:length(U2.Legend)
                if m == length(U2.Legend)
                    U2 = lab_plot_legend(U2.Legend{m},U2,1);
                else
                    U2 = lab_plot_legend(U2.Legend{m},U2);
                end
            end
            Lflag = true;
        end
        set(gcf,'UserData',U2);
    end
     
    function show_legend
        for j = 1:length(PLOT.Legend)
            Table{j,1} = lab_color2html(PLOT.Legend(j).Color,PLOT.Legend(j).Text); %#ok<AGROW>
        end
        lab_table_dialog(Table,{'Name'},'Legend',-1,[],[1 1 1]);
    end

     function adjust_threshold
         value = get(Hslider,'Value');
         numel = length(settings.handleE);
         val = round(value * numel);
         if val > 0
             set(settings.handleE(1:val),'visible','off');
         end
         if val < numel
             set(settings.handleE(val+1:numel),'visible','on');
         end
     end
     
     function adjust_threshold_menu
         Tset.pvalue = get(Hslider,'Value') * 100;
         Tset.value = round(1000 * ((Tset.pvalue/100) * (MatrixMax - MatrixMin) + MatrixMin)) / 1000;
         
         numel = length(settings.handleE);
         Prompt = {'Threshold','value',''};
         Formats.type = 'edit';
         Formats.format = 'float';
         Formats.size = 40;
         Formats.limits = [MatrixMin MatrixMax];
         Formats.callback = {@set_pvalue,'@ALL','@ALL',MatrixMin,MatrixMax};
         
         Prompt(end+1,:) = {' = ','pvalue','%'};
         Formats(end+1,1).type = 'edit';
         Formats(end,1).format = 'float';
         Formats(end,1).size = 40;
         Formats(end,1).limits = [0 100];
         Formats(end,1).callback = {@set_value,'@ALL','@ALL',MatrixMin,MatrixMax};  
         
         Tset = inputsdlg(Prompt,'Set Threshold',Formats,Tset,2);
         val = round((Tset.pvalue/100) * numel);
         if val > 0
             set(settings.handleE(1:val),'visible','off');
         end
         if val < numel
             set(settings.handleE(val+1:numel),'visible','on');
         end
         set(Hslider,'Value',Tset.pvalue/100);
     end
end

function Tset = set_value(Tset,Mmin,Mmax)
    if ~isempty(Tset.pvalue)
        Tset.value = round(1000 * ((Tset.pvalue/100) * (Mmax - Mmin) + Mmin)) / 1000;
    else
        Tset.value = [];
    end
end

function Tset = set_pvalue(Tset,Mmin,Mmax)
    if ~isempty(Tset.value)
        Tset.pvalue = round(1000 * ((Tset.value - Mmin) / (Mmax - Mmin)) * 100) / 1000;
    else
        Tset.pvalue = [];
    end
end