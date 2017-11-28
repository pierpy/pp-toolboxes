function lab_plot_IS_plotsingle(PlotFaces,PlotVertices,facecolor,PlotFacesR,PlotVerticesR,facecolorR,PlotFacesL,PlotVerticesL,facecolorL,Matrix,Valpha,PLOT,cfg)

Modes = fieldnames(PLOT);
Ptitle = '';
Legend = {};
for M = 1:length(Modes)
    if isfield(PLOT.(Modes{M}),'Name')
        if isempty(Ptitle)
            Ptitle =  PLOT.(Modes{M}).Name;
        else
            Ptitle =  [Ptitle PLOT.(Modes{M}).Name]; %#ok<AGROW>
        end
    end
    if isfield(PLOT.(Modes{M}),'Legend') & isfield(PLOT.(Modes{M}).Legend,'Mode')
        Legend{end+1,1} = PLOT.(Modes{M}).Legend; %#ok<AGROW>
    end
end
cfg.backgroundcolor = [1 1 1];
cfg.textcolor = [0 0 0];
figure1 = figure('Color',cfg.backgroundcolor,'InvertHardCopy','off','NumberTitle','off','Menubar','none');
%figure1 = figure('Color',cfg.backgroundcolor,'InvertHardCopy','off','NumberTitle','off','Renderer','zbuffer');
[~,~,~,DATA_fileS] = lab_filename(cfg.DATA_file);
if ~isempty(Ptitle)
    set(gcf,'Name',[DATA_fileS '_' Ptitle]);
else
    set(gcf,'Name',DATA_fileS);
end
Thandle = [];
Hlight = [];
Lflag = false;
if isfield(cfg,'nosingle') & isfield(cfg,'DATA_file') & ~isempty(cfg.DATA_file) & cfg.Store == true
    dosplit = true;
    show_brain;
    close(figure1);
else
    dosplit = false;
    show_brain;
end

% extra scripts
    function change_view
        if dosplit == true
            dosplit = false;
        else
            dosplit = true;
        end
        U2 = get(gcf,'UserData');
        U2 = lab_plot_legend('off',U2);
        set(gcf,'UserData',U2);
        Lflag = false;
        show_brain
    end

    function show_brain
        clf(figure1);
        Thandle = [];
        Ahandle = [];
        m1 = uimenu(figure1,'Label','File');
        uimenu(m1,'Label','Save Picture','Callback',@(~,~)lab_print_figure);
        uimenu(m1,'Label','Copy Picture','Callback','print -dmeta;');
        uimenu(m1,'Label','Close','Callback','close;');
        m2 = uimenu(figure1,'Label','Edit');
        uimenu(m2,'Label','Black/White Background','Callback',@(~,~)set_background);
        uimenu(m2,'Label','Change view','Callback',@(~,~)change_view);
        uimenu(m2,'Label','Rotate on/off','Callback',@(~,~)rotate3d);
        m3 = uimenu(figure1,'Label','View');
        uimenu(m3,'Label','Top','Callback',@(~,~)view_top);
        uimenu(m3,'Label','Bottom','Callback',@(~,~)view_bottom);
        uimenu(m3,'Label','Left','Callback',@(~,~)view_left);
        uimenu(m3,'Label','Right','Callback',@(~,~)view_right);
        uimenu(m3,'Label','Front','Callback',@(~,~)view_front);
        uimenu(m3,'Label','Back','Callback',@(~,~)view_back);
        if ~isempty(Legend)
            uimenu(m3,'Label','Toggle Legend','Callback',@(~,~)toggle_legend,'Separator','on');
            uimenu(m3,'Label','Edit Legend','Callback',@(~,~)toggle_legend2);
        end
        if isfield(cfg.PLOT,'Channels')
            uimenu(m3,'Label','Channels-Names','Callback',@(~,~)show_legend);
        end
        if dosplit == true
            if isempty(Matrix)
                Ahandle(end+1) = lab_show_brain(figure1,PlotFaces,PlotVertices,facecolor,Matrix,Valpha,'L',[2 3 1],0);
                Ahandle(end+1) = lab_show_brain(figure1,PlotFaces,PlotVertices,facecolor,Matrix,Valpha,'T',[2 3 2],0);
                Ahandle(end+1) = lab_show_brain(figure1,PlotFaces,PlotVertices,facecolor,Matrix,Valpha,'R',[2 3 3],0);
                Ahandle(end+1) = lab_show_brain(figure1,PlotFacesL,PlotVerticesL,facecolorL,Matrix,Valpha,'R',[2 3 4],0);
                Ahandle(end+1) = lab_show_brain(figure1,PlotFaces,PlotVertices,facecolor,Matrix,Valpha,'B',[2 3 5],0);
                Ahandle(end+1) = lab_show_brain(figure1,PlotFacesR,PlotVerticesR,facecolorR,Matrix,Valpha,'L',[2 3 6],0);
            else
                if unique(Matrix.edges) == 2
                    threshold1 = 0;
                    threshold2 = [];
                    threshold3 = [];
                elseif length(Matrix.edges) <= 230
                    threshold1 = 0;
                    threshold2 = 0.5;
                    threshold3 = 0.75;
                else
                    threshold1 = 0.75;
                    threshold2 = 0.875;
                    threshold3 = 0.9688;
                end
                if ~isempty(threshold3)
                    Ahandle(end+1) = lab_show_brain(figure1,PlotFaces,PlotVertices,facecolor,Matrix,Valpha,'T',[2 3 1],0,threshold1);
                    Ahandle(end+1) = lab_show_brain(figure1,PlotFaces,PlotVertices,facecolor,Matrix,Valpha,'T',[2 3 2],0,threshold2);
                    Ahandle(end+1) = lab_show_brain(figure1,PlotFaces,PlotVertices,facecolor,Matrix,Valpha,'T',[2 3 3],0,threshold3);
                    Ahandle(end+1) = lab_show_brain(figure1,PlotFaces,PlotVertices,facecolor,Matrix,Valpha,'L',[2 3 4],0,threshold1);
                    Ahandle(end+1) = lab_show_brain(figure1,PlotFaces,PlotVertices,facecolor,Matrix,Valpha,'L',[2 3 5],0,threshold2);
                    Ahandle(end+1) = lab_show_brain(figure1,PlotFaces,PlotVertices,facecolor,Matrix,Valpha,'L',[2 3 6],0,threshold3);
                    Thandle = uicontrol('Style', 'text','String',['min percent: ' num2str(threshold1*100) '%  ' num2str(threshold2*100) '%  ' num2str(threshold3*100) '%'], ...
                        'Units','normalized','Position',[0 0 0.25 0.05],'BackgroundColor',cfg.backgroundcolor, ...
                        'ForegroundColor',cfg.textcolor);
                elseif ~isempty(threshold2)
                    Ahandle(end+1) = lab_show_brain(figure1,PlotFaces,PlotVertices,facecolor,Matrix,Valpha,'T',[2 3 1],0,threshold1);
                    Ahandle(end+1) = lab_show_brain(figure1,PlotFaces,PlotVertices,facecolor,Matrix,Valpha,'L',[2 3 2],0,threshold1);
                    Ahandle(end+1) = lab_show_brain(figure1,PlotFaces,PlotVertices,facecolor,Matrix,Valpha,'F',[2 3 3],0,threshold1);
                    Ahandle(end+1) = lab_show_brain(figure1,PlotFaces,PlotVertices,facecolor,Matrix,Valpha,'T',[2 3 4],0,threshold2);
                    Ahandle(end+1) = lab_show_brain(figure1,PlotFaces,PlotVertices,facecolor,Matrix,Valpha,'L',[2 3 5],0,threshold2);
                    Ahandle(end+1) = lab_show_brain(figure1,PlotFaces,PlotVertices,facecolor,Matrix,Valpha,'F',[2 3 6],0,threshold2);
                    Thandle = uicontrol('Style', 'text','String',['min percent: ' num2str(threshold1*100) '%  ' num2str(threshold2*100) '%'], ...
                        'Units','normalized','Position',[0 0 0.2 0.05],'BackgroundColor',cfg.backgroundcolor, ...
                        'ForegroundColor',cfg.textcolor);
                else
                    Ahandle(end+1) = lab_show_brain(figure1,PlotFaces,PlotVertices,facecolor,Matrix,Valpha,'T',[1 3 1],0,threshold1);
                    Ahandle(end+1) = lab_show_brain(figure1,PlotFaces,PlotVertices,facecolor,Matrix,Valpha,'L',[1 3 2],0,threshold1);
                    Ahandle(end+1) = lab_show_brain(figure1,PlotFaces,PlotVertices,facecolor,Matrix,Valpha,'F',[1 3 3],0,threshold1);
                    Thandle = uicontrol('Style', 'text','String',['min precent: ' num2str(threshold1*100) '%'], ...
                        'Units','normalized','Position',[0 0 0.2 0.05],'BackgroundColor',cfg.backgroundcolor, ...
                        'ForegroundColor',cfg.textcolor);
                end
            end
            
            % plot legend
            if ~isempty(Legend)
                U = get(gcf,'UserData');
                U.Ahandle = Ahandle;
                U.Fhandle = figure1;
                Lflag = false;
                set(gcf,'UserData',U);
            end
            
            % save picture
            if isfield(cfg,'DATA_file') & cfg.Store == true
                lab_print_figure(cfg.DATA_file,figure1);
            end
        else
            if isempty(Matrix)
                [Ahandle,Hlight] = lab_show_brain(figure1,PlotFaces,PlotVertices,facecolor,Matrix,Valpha,'T',1,0);
                set(Ahandle,'FontSize',8);
                Hlight(2) = light('Position',[0 0 -1]);
                rotate3d on
            else
                [Ahandle,Hlight] = lab_show_brain(figure1,PlotFaces,PlotVertices,facecolor,Matrix,Valpha,'T',1,1);
                Hlight(2) = light('Position',[0 0 -1]);
                rotate3d on
            end
            
            % plot legend
            if ~isempty(Legend)
                U = get(gcf,'UserData');
                U.Ahandle = Ahandle;
                U.Fhandle = figure1;
                for L = 1:length(Legend)
                    U = lab_plot_legend(Legend{L},U);
                end
                Lflag = true;
                set(gcf,'UserData',U);
            end
            
            % save picture
            if isfield(cfg,'DATA_file') & cfg.Store == true
                lab_print_figure(cfg.DATA_file,figure1);
            end
        end
        
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
        if ~isempty(Thandle)
            for T = 1:length(Thandle)
                tmp = get(Thandle(T));
                if isfield(tmp,'ForegroundColor')
                    set(Thandle(T),'ForegroundColor',textcolor);
                end
                if isfield(tmp,'XColor')
                    set(Thandle(T),'XColor',textcolor);
                end
            end
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
    function toggle_legend
        U2 = get(gcf,'UserData');
        if Lflag == true
            U2 = lab_plot_legend('off',U2);
            Lflag = false;
        else
            for m = 1:length(Legend)
                U2 = lab_plot_legend(Legend{m},U2);
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
            for m = 1:length(Legend)
                if m == length(Legend)
                    U2 = lab_plot_legend(Legend{m},U2,1);
                else
                    U2 = lab_plot_legend(Legend{m},U2);
                end
            end
            Lflag = true;
        end
        set(gcf,'UserData',U2);
    end
    function show_legend
        for j = 1:length(cfg.PLOT)
            Table{j,1} = lab_color2html(cfg.PLOT(j).Color1,cfg.PLOT(j).Name); %#ok<AGROW>
        end
        lab_table_dialog(Table,{'Name'},'Legend',-1,[],[1 1 1]);
    end
end