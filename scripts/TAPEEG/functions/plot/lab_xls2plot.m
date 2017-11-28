% Plot boxplot / scatter / q-q-plot
%
% settings2 = lab_xls2plot(settings)
%
%     Input         xls-file with results/outcomes (see lab_read_statistics)
%
% Written by F. Hatz 2013

function settings2 = lab_xls2plot(settings)

settings2 = [];
if ~exist('settings','var')
    settings = [];
end
[data,header,result,~,settings] = lab_read_statistics(settings,1,0,1,0,1);
if isempty(data)
    return
end
filenameS = settings.file;
filepath = settings.path;
if isempty(data)
    return
end
varnames = {};
for V = 1:size(data,2)
    tmpV = unique(data(:,V),'stable');
    varnames{end+1,1} = cellstr(num2str(tmpV(:))); %#ok<AGROW>
end
if isfield(settings,'ResultNames')
    varnames = cat(1,varnames,settings.ResultNames(:));
elseif ~isempty(result)
    for V = 1:size(result,2)
        tmpV = unique(result(:,V),'stable');
        varnames{end+1,1} = cellstr(num2str(tmpV(:))); %#ok<AGROW>
    end
end
vars = [header.vars(:);header.result(:)];
data = cat(1,data',result');
varvalid = 1:size(data,2);

do_settings;
fig = figure;
do_figure;
if isempty(settings)
    close;
    return
end
filesave = plot_xls;
lab_print_figure(fullfile(filepath,filesave),fig,settings.fileformat);

    function do_figure
        set(gcf,'Color',[1 1 1],'Menubar','none','Toolbar','none');
        m1 = uimenu(gcf,'Label','File');
        uimenu(m1,'Label','Save Picture','Callback',@(~,~)lab_print_figure);
        uimenu(m1,'Label','Copy Picture','Callback','print -dmeta;');
        uimenu(m1,'Label','Close','Callback','close;');
        m2 = uimenu(gcf,'Label','Plot');
        uimenu(m2,'Label','Replace','Callback',@(~,~)do_plot_replace);
        uimenu(m2,'Label','New','Callback',@(~,~)do_plot_new);
        m2 = uimenu(gcf,'Label','Edit');
        uimenu(m2,'Label','Color/Gray','Callback',@(~,~)switch_colormode);
        uimenu(m2,'Label','Y-Range','Callback',@(~,~)Ylimits);
        uimenu(m2,'Label','X-Range','Callback',@(~,~)Xlimits);
    end

    function do_plot_replace(skipset)
        if ~exist('skipset','var') | skipset ~= 1
            do_settings;
        end
        if ~isempty(settings)
            clf;
            do_figure;
            filesave = plot_xls;
            lab_print_figure(fullfile(filepath,filesave),fig,settings.fileformat);
        end
    end

    function do_plot_new
        fig = figure;
        do_figure;
        do_settings;
        if isempty(settings)
            close;
            return
        end
        filesave = plot_xls;
        lab_print_figure(fullfile(filepath,filesave),fig,settings.fileformat);
    end

    function do_settings
        [settings,skipprocessing] = lab_set_xls2plot(settings,data,vars,varnames,[filenameS '.xlsx']);
        if skipprocessing == 1
            settings = [];
            return
        else
            settings2 = settings;
        end
    end

    function filenameout = plot_xls
        % find valid trials
        groupvalid = 1:size(data,2);
        if ~isempty(settings.xvar)
            groupvalid = setdiff(groupvalid,find(isnan(data(settings.xvar,:))));
            groupvalid = setdiff(groupvalid,find(isempty(data(settings.xvar,:))));
        end
        varvalid = 1:size(data,2);
        for i = 1:length(settings.yvar)
            varvalid = setdiff(varvalid,find(isnan(data(settings.yvar(i),:))));
            varvalid = setdiff(varvalid,find(isempty(data(settings.yvar(i),:))));
        end
        varvalid = intersect(varvalid,groupvalid);
        if ~isempty(settings.xvar)
            tmp = data(settings.xvar,varvalid);
            [~,tmp] = sort(tmp);
            varvalid = varvalid(tmp);
        end
                
        % do plot
        if strcmp(settings.mode,'Boxplot')
            filenameout = do_boxplot;
        elseif strcmp(settings.mode,'Q-Q-Plot')
            filenameout = do_q_q_plot;
        elseif strcmp(settings.mode,'Q-Q-Plot-Normal')
            filenameout = do_q_q_plot_normal;
        elseif strcmp(settings.mode,'Scatter')
            filenameout = do_scatter;
        end
    end

    function filenameout = do_boxplot
        plotdata = [];
        plotgroup = [];
        p = zeros(1,length(settings.yvar));
        p2 = [];
        p3 = [];
        for i = 1:length(settings.yvar)
            if ~isempty(settings.xvar)
                plotdata = [plotdata data(settings.yvar(i),varvalid)]; %#ok<AGROW>
                plotgroup = [plotgroup cat(1,ones(1,size(data(settings.xvar,varvalid),2))*i,data(settings.xvar,varvalid))]; %#ok<AGROW>
                Rvars = unique(data(settings.xvar,varvalid));
                Rvars = sort(Rvars);
                if length(Rvars) == 2
                    for j = 1:2
                        dat{j} = data(settings.yvar(i),data(settings.xvar,:)==Rvars(j)); %#ok<AGROW>
                        if settings.ylog(i) == true
                            dat{j} = log(dat{j}); %#ok<AGROW>
                        end
                    end
                    if strcmp(settings.statistics,'T-test')
                        [~,p(i)] = ttest2(dat{1}',dat{2}');
                    elseif strcmp(settings.statistics,'Mann-Whitney-U')
                        p(i) = ranksum(dat{1}',dat{2}');
                    end
                    clearvars dat
                else
                    datatest{1} = data(settings.yvar(i),varvalid)';
                    if settings.ylog(i) == true
                        datatest{1} = log(datatest{1});
                    end
                    datatest{2} = data(settings.xvar,varvalid)';
                    if strcmp(settings.statistics,'Anova')
                        p(i) = anova1(datatest{1},datatest{2},'off');
                    elseif strcmp(settings.statistics,'KruskalWallis')
                        p(i) = kruskalwallis(datatest{1},datatest{2},'off');
                    end
                    for j = 2:length(Rvars)
                        dat{1} = data(settings.yvar(i),data(settings.xvar,:)==Rvars(j-1));
                        dat{2} = data(settings.yvar(i),data(settings.xvar,:)==Rvars(j));
                        if settings.ylog(i) == true
                            dat{1} = log(dat{1});
                            dat{2} = log(dat{2});
                        end
                        if strcmp(settings.statistics,'Anova')
                            [~,p2(i,j-1)] = ttest2(dat{1}',dat{2}'); %#ok<AGROW>
                        elseif strcmp(settings.statistics,'KruskalWallis')
                            p2(i,j-1) = ranksum(dat{1}',dat{2}'); %#ok<AGROW>
                        end
                        clearvars dat
                    end
                    if length(Rvars) == 3
                        dat{1} = data(settings.yvar(i),data(settings.xvar,:)==Rvars(1));
                        dat{2} = data(settings.yvar(i),data(settings.xvar,:)==Rvars(3));
                        if settings.ylog(i) == true
                            dat{1} = log(dat{1});
                            dat{2} = log(dat{2});
                        end
                        if strcmp(settings.statistics,'Anova')
                            [~,p3(i,1)] = ttest2(dat{1}',dat{2}'); %#ok<AGROW>
                        elseif strcmp(settings.statistics,'KruskalWallis')
                            p3(i,1) = ranksum(dat{1}',dat{2}'); %#ok<AGROW>
                        end
                        clearvars dat
                    end
                end
            else
                plotdata = [plotdata data(settings.yvar(i),varvalid)']; %#ok<AGROW>
            end
        end
        if isempty(settings.xvar)
            labels = settings.yname';
        elseif length(settings.yvar) > 1
            labels = repmat(settings.xname',1,length(settings.yvar));
        else
            labels = settings.xname';
        end
        warning off %#ok<WNOFF>
        if ~isempty(plotgroup)
            C = boxplot(plotdata',plotgroup','labels',labels,'Colors','k','Symbol','+k','factorgap',[8 0]);
        else
            C = boxplot(plotdata,'labels',labels,'Colors','k','Symbol','+k');
        end
        for i = 1:size(C,2)
            tmp = get(C(1,i),'XData');
            Cx(i) = tmp(1); %#ok<AGROW>
        end
        clearvars tmp
        if ~isempty(settings.xvar)
            for i = 1:length(settings.yvar)
                Cxvar(i) = mean(Cx((i-1)*length(Rvars)+1:i*length(Rvars))); %#ok<AGROW>
            end
        end
        if settings.ylog(i) == true
            set(gca,'YScale','log');
        end
        if ~isempty(settings.Limits) & isfield(settings.Limits,'y1') & ~isempty(settings.Limits.y1) & ...
                isfield(settings.Limits,'y2') & ~isempty(settings.Limits.y2)
            set(gca,'YLim',[settings.Limits.y1 settings.Limits.y2]);
        end
        
        text_h = findobj(gca,'Type','text');
        if length(labels) > 8
            for cnt = 1:length(labels)
                set(text_h(end-cnt+1),'Rotation',45, ...
                    'String',labels{cnt},'HorizontalAlignment','right');
            end
        end
        mt = uimenu(gcf,'Label','RotateText');
        uimenu(mt,'Label','0 grad','Callback',@(~,~)set_X_Angle(text_h,labels,0));
        uimenu(mt,'Label','15 grad','Callback',@(~,~)set_X_Angle(text_h,labels,15));
        uimenu(mt,'Label','30 grad','Callback',@(~,~)set_X_Angle(text_h,labels,30));
        uimenu(mt,'Label','45 grad','Callback',@(~,~)set_X_Angle(text_h,labels,45));
        uimenu(mt,'Label','60 grad','Callback',@(~,~)set_X_Angle(text_h,labels,60));
        uimenu(mt,'Label','90 grad','Callback',@(~,~)set_X_Angle(text_h,labels,90));
        if length(settings.yvar) > 1
            if ~isempty(settings.xvar)
                T = get(gca);
                if settings.ylog(end) == true
                    Ycorr = (T.YLim(2)-T.YLim(1));
                else
                    Ycorr = (T.YLim(2)-T.YLim(1))*0.03;
                end
                T.YLim(2) = T.YLim(2) + 5*Ycorr;
                set(gca,'YLim',T.YLim);
                Yv = T.YLim(2) - Ycorr;
                Yp = T.YLim(2) - 2.6*Ycorr;
                Yp2 = T.YLim(2) - 4.2*Ycorr;
                if ~isempty(p3)
                    Yp2 = T.YLim(2) - 5.2*Ycorr;
                    Yp3 = T.YLim(2) - 4*Ycorr;
                end
                for i = 1:length(settings.yvar)
                    if ~isempty(p3)
                        C1 = Cx((i-1)*length(Rvars)+1);
                        C2 = Cx((i-1)*length(Rvars)+3);
                        hold on
                        plot([C1+0.01 C1+0.01],[Yp3-Ycorr/2 Yp3+Ycorr/2],'Color',[0 0 0]);
                        plot([C2-0.01 C2-0.01],[Yp3-Ycorr/2 Yp3+Ycorr/2],'Color',[0 0 0]);
                        plot([C1+0.01 C2-0.01],[Yp3 Yp3],'Color',[0 0 0]);
                        Ctmp = mean([C1 C2]);
                        if p3(i,1) >= 0.01
                            text(Ctmp,Yp3,num2str(p3(i,1),'%1.2f'),'HorizontalAlignment','center','FontSize',8,'BackgroundColor',[1 1 1]);
                        elseif p3(i,1) > 0
                            text(Ctmp,Yp3,'<0.01','HorizontalAlignment','center','FontSize',8,'BackgroundColor',[1 1 1]);
                        end
                    end
                    for j = 1:size(p2,2)
                        C1 = Cx((i-1)*length(Rvars)+j);
                        C2 = Cx((i-1)*length(Rvars)+j+1);
                        hold on
                        plot([C1+0.01 C1+0.01],[Yp2-Ycorr/2 Yp2+Ycorr/2],'Color',[0 0 0]);
                        plot([C2-0.01 C2-0.01],[Yp2-Ycorr/2 Yp2+Ycorr/2],'Color',[0 0 0]);
                        plot([C1+0.01 C2-0.01],[Yp2 Yp2],'Color',[0 0 0]);
                        Ctmp = mean([C1 C2]);
                        if p2(i,j) >= 0.01
                            text(Ctmp,Yp2,num2str(p2(i,j),'%1.2f'),'HorizontalAlignment','center','FontSize',8,'BackgroundColor',[1 1 1]);
                        elseif p2(i,j) > 0
                            text(Ctmp,Yp2,'<0.01','HorizontalAlignment','center','FontSize',8,'BackgroundColor',[1 1 1]);
                        end
                    end
                    text(Cxvar(i),Yv,settings.yname{i},'HorizontalAlignment','center');
                    if p(i) >= 0.01
                        text(Cxvar(i),Yp,['p = ' num2str(p(i),'%1.2f')],'HorizontalAlignment','center','FontSize',8);
                    elseif p(i) > 0
                        text(Cxvar(i),Yp,'p < 0.01','HorizontalAlignment','center','FontSize',8);
                    end
                end
            end
            if ~isempty(settings.title)
                title(settings.title);
            end
        else
            if ~isempty(settings.xvar)
                T = get(gca);
                if settings.ylog(end) == true
                    Ycorr = (T.YLim(2)-T.YLim(1));
                else
                    Ycorr = (T.YLim(2)-T.YLim(1))*0.03;
                end
                T.YLim(2) = T.YLim(2) + 3.8*Ycorr;
                set(gca,'YLim',T.YLim);
                Yp2 = T.YLim(2) - 2*Ycorr;
                if ~isempty(p3)
                    Yp2 = T.YLim(2) - 3.4*Ycorr;
                    Yp3 = T.YLim(2) - 2*Ycorr;
                end
                if ~isempty(p3)
                    C1 = Cx((i-1)*length(Rvars)+1);
                    C2 = Cx((i-1)*length(Rvars)+3);
                    hold on
                    plot([C1+0.01 C1+0.01],[Yp3-Ycorr/2 Yp3+Ycorr/2],'Color',[0 0 0]);
                    plot([C2-0.01 C2-0.01],[Yp3-Ycorr/2 Yp3+Ycorr/2],'Color',[0 0 0]);
                    plot([C1+0.01 C2-0.01],[Yp3 Yp3],'Color',[0 0 0]);
                    Ctmp = mean([C1 C2]);
                    if p3(i,1) >= 0.01
                        text(Ctmp,Yp3,num2str(p3(i,1),'%1.2f'),'HorizontalAlignment','center','FontSize',8,'BackgroundColor',[1 1 1]);
                    elseif p3(i,1) > 0
                        text(Ctmp,Yp3,'<0.01','HorizontalAlignment','center','FontSize',8,'BackgroundColor',[1 1 1]);
                    end
                end
                for j = 1:size(p2,2)
                    C1 = Cx((i-1)*length(Rvars)+j);
                    C2 = Cx((i-1)*length(Rvars)+j+1);
                    hold on
                    plot([C1+0.01 C1+0.01],[Yp2-Ycorr/2 Yp2+Ycorr/2],'Color',[0 0 0]);
                    plot([C2-0.01 C2-0.01],[Yp2-Ycorr/2 Yp2+Ycorr/2],'Color',[0 0 0]);
                    plot([C1+0.01 C2-0.01],[Yp2 Yp2],'Color',[0 0 0]);
                    Ctmp = mean([C1 C2]);
                    if p2(i,j) >= 0.01
                        text(Ctmp,Yp2,num2str(p2(i,j),'%1.2f'),'HorizontalAlignment','center','FontSize',8,'BackgroundColor',[1 1 1]);
                    elseif p2(i,j) > 0
                        text(Ctmp,Yp2,'<0.01','HorizontalAlignment','center','FontSize',8,'BackgroundColor',[1 1 1]);
                    end
                end
            end
            
            if p(1) >= 0.001
                plottitle = [settings.yname{1} ' (p = ' num2str(p(1),'%1.3f') ')'];
            elseif p(1) > 0
                plottitle = [settings.yname{1} ' (p < 0.001)'];
            else
                plottitle = settings.yname{1};
            end
            if ~isempty(settings.title)
                title([plottitle ' - ' settings.title]);
            else
                title(plottitle);
            end
        end
        varnameout = settings.yname{1};
        if length(settings.yname) <= 4
            for i = 2:length(settings.yname)
                varnameout = [varnameout '_' settings.yname{i}]; %#ok<AGROW>
            end
        else
            varnameout = [settings.yname{1} '_to_' settings.yname{end}];
        end
        set(gcf,'Name',[varnameout ' - ' vars{settings.xvar}],'NumberTitle','off');
        filenameout = [filenameS '_' varnameout '_' vars{settings.xvar} '_boxplot.jpg'];
        filenameout = regexprep(filenameout,{':',',',';','\','/'},'');
        warning on %#ok<WNON>
    end

    function filenameout = do_q_q_plot
        qqplot(data(settings.xvar,varvalid)',data(settings.yvar(1),varvalid)');
        H = kstest2(data(settings.xvar,varvalid)',data(settings.yvar(1),varvalid)');
        Title = ['Q-Q-Plot ' regexprep(settings.xname{1},'_',' ') ' versus ' ...
            regexprep(settings.yname{1},'_',' ') ' (KS: ' num2str(H) ')'];
        title(Title);
        set(gcf,'Name',Title,'NumberTitle','off');
        xlabel(settings.xname{1});
        ylabel(settings.yname{1});
        filenameout = [filenameS '_' settings.xname{1} '_' settings.yname{1} '_qqplot.jpg'];
        filenameout = regexprep(filenameout,{':',',',';','\'},'');
        if strcmp(settings.colormode,'gray')
            tmp = get(gca,'Children');
            set(tmp(1),'MarkerEdgeColor',[0.5 0.5 0.5]);
            set(tmp(2),'Color',[0 0 0]);
            set(tmp(3),'Color',[0 0 0]);
        end
    end

    function filenameout = do_q_q_plot_normal
        qqplot(data(settings.xvar,varvalid)');
        H = kstest(data(settings.xvar,varvalid)');
        Title = ['Q-Q-Plot ' regexprep(settings.xname{1},'_',' ') ' (KS: ' num2str(H) ')'];
        title(Title);
        set(gcf,'Name',Title,'NumberTitle','off');
        xlabel(settings.xname{1});
        ylabel('Normal distribution');
        filenameout = [filenameS '_' settings.xname{1} '_qqplot.jpg'];
        filenameout = regexprep(filenameout,{':',',',';','\'},'');
        if strcmp(settings.colormode,'gray')
            tmp = get(gca,'Children');
            set(tmp(1),'MarkerEdgeColor',[0.5 0.5 0.5]);
            set(tmp(2),'Color',[0 0 0]);
            set(tmp(3),'Color',[0 0 0]);
        end
    end

    function filenameout = do_scatter
        p = zeros(1,length(settings.yvar));
        if strcmp(settings.statistics,'Kendall') | strcmp(settings.statistics,'Spearman') | strcmp(settings.statistics,'Pearson')
            for i = 1:length(settings.yvar)
                if settings.ylog(i) == true
                    datatmp = log(data(settings.yvar(i),varvalid)');
                else
                    datatmp = data(settings.yvar(i),varvalid)';
                end
                if settings.xlog == true
                    grouptmp = log(data(settings.xvar,varvalid)');
                else
                    grouptmp = data(settings.xvar,varvalid)';
                end
                [~,p(i)]=corr(datatmp,grouptmp,'type',settings.statistics);
            end
        end
        varnames = {};
        for i = 1:length(settings.yvar)
            if p(i) >= 0.001
                varnames = [varnames cellstr([settings.yname{i} ' (p = ' num2str(p(i),'%1.3f') ')'])]; %#ok<AGROW>
            elseif p(i) > 0
                varnames = [varnames cellstr([settings.yname{i} ' (p < 0.001)'])]; %#ok<AGROW>
            else
                varnames = [varnames cellstr(settings.yname{i})]; %#ok<AGROW>
            end
        end
        if strcmp(settings.colormode,'color')
            tmp = hsv(length(settings.yvar));
            tmp = tmp(randperm(length(settings.yvar)),:);
            colormap = cat(1,[0 0 0],tmp);
        else
            tmp = gray(round((length(settings.yvar)+1)));
            tmp(2:end,:) = tmp(randperm(length(settings.yvar))+1,:);
            colormap = tmp;
        end
        plotgroup = [];
        plotdata = [];
        if min(settings.ylog) == 1
            varlogall = true;
        else
            varlogall = false;
        end
        for i = 1:length(settings.yvar)
            if settings.ylog(i) == true & varlogall == false
                datatmp = log(data(settings.yvar(i),varvalid)');
            else
                datatmp = data(settings.yvar(i),varvalid)';
            end
            if isempty(settings.Marker)
                settings.Marker.size = 36;
                settings.Marker.marker = 'o';
            end
            scatter(data(settings.xvar,varvalid)',datatmp,settings.Marker.size,colormap(i,:), ...
                settings.Marker.marker,'DisplayName',settings.yname{i});
            hold on
            plotgroup = cat(1,plotgroup,data(settings.xvar,varvalid)');
            plotdata = cat(1,plotdata,data(settings.yvar(i),varvalid)');
        end
        if ~isempty(settings.Limits) & isfield(settings.Limits,'y1') & ~isempty(settings.Limits.y1) & ...
                isfield(settings.Limits,'y2') & ~isempty(settings.Limits.y2)
            set(gca,'YLim',[settings.Limits.y1 settings.Limits.y2]);
        end
        if ~isempty(settings.Limits) & isfield(settings.Limits,'x1') & ~isempty(settings.Limits.x1) & ...
                isfield(settings.Limits,'x2') & ~isempty(settings.Limits.x2)
            set(gca,'XLim',[settings.Limits.x1 settings.Limits.x2]);
        end
        if varlogall == true
            set(gca,'YScale','log');
        end
        if settings.xlog == true
            set(gca,'XScale','log');
        end
        if settings.XequalY == true
            XLim = get(gca,'Xlim');
            YLim = get(gca,'YLim');
            XLim(1) = max(XLim(1),YLim(1));
            XLim(2) = max(XLim(2),YLim(2));
            set(gca,'XLim',XLim);
            set(gca,'YLim',XLim);
        end
        pfit = polyfit(plotgroup,plotdata,1);
        x(1) = min(plotgroup);
        x(2) = max(plotgroup);
        y = polyval(pfit,x);
        hold on
        if length(settings.yvar) > 1 & strcmp(settings.colormode,'color')
            plot(x,y,'r-');
        else
            plot(x,y,'k-');
        end
        xlabel(settings.xname{1});
        if length(settings.yvar) == 1
            ylabel(settings.yname{1});
            if strcmp(settings.statistics,'None')
                if ~isempty(settings.title)
                    title(settings.title);
                end
            else
                if p(1) < 0.001
                    plottitle = 'p < 0.001 ';
                else
                    plottitle = ['p = ' num2str(p(1),'%1.3f')];
                end
                if ~isempty(settings.title)
                    title([settings.title ' (' plottitle ')']);
                else
                    title(plottitle);
                end
            end
        else
            if ~isempty(settings.title)
                title(settings.title);
            end
            if isempty(settings.Legend)
                settings.Legend.border = true;
                settings.Legend.place = 'NorthEast';
            end
            legend(varnames,'Location',settings.Legend.place);
            if settings.Legend.border == false
                legend('boxoff');
            end
        end
        set(gcf,'Name',[settings.title ' - ' settings.xname{1}],'NumberTitle','off');
        tmp = settings.xname{1};
        for i = 1:length(settings.yname)
            tmp = [tmp '_' settings.yname{i}]; %#ok<AGROW>
        end
        filenameout = [filenameS '_' settings.title '_' tmp '_scatter.jpg'];
        filenameout = regexprep(filenameout,{':',',',';','\'},'');
    end

    function switch_colormode
        if strcmp(settings.colormode,'color')
            settings.colormode = 'gray';
        else
            settings.colormode = 'color';
        end
        do_plot_replace(1);
    end

    function Xlimits
        if strcmp(settings.mode,'Scatter')
            XLim = get(gca,'XLim');
            settings.Limits.x1 = XLim(1);
            settings.Limits.x2 = XLim(2);
            Prompt = {'Min value','x1';'Max value','x2'};
            Formats.type = 'edit';
            Formats.format = 'float';
            Formats.size = 30;
            Formats.limits = [-inf inf];
            Formats = [Formats Formats];
            [settings.Limits,Cancelled] = inputsdlg(Prompt,'X-Range',Formats,settings.Limits);
            if Cancelled ~= 1
                XLim = [settings.Limits.x1 settings.Limits.x2];
                set(gca,'XLim',XLim);
            end
            clearvars Formats Prompt
        end
    end

    function Ylimits
        if strcmp(settings.mode,'Boxplot') | strcmp(settings.mode,'Scatter')
            YLim = get(gca,'YLim');
            settings.Limits.y1 = YLim(1);
            settings.Limits.y2 = YLim(2);
            Prompt = {'Min value','y1';'Max value','y2'};
            Formats.type = 'edit';
            Formats.format = 'float';
            Formats.size = 30;
            Formats.limits = [-inf inf];
            Formats = [Formats Formats];
            [settings.Limits,Cancelled] = inputsdlg(Prompt,'Y-Range',Formats,settings.Limits);
            if Cancelled ~= 1
                YLim = [settings.Limits.y1 settings.Limits.y2];
                set(gca,'YLim',YLim);
            end
            clearvars Formats Prompt
        end
    end
end

function set_X_Angle(text_h,labels,angle)
    if length(text_h) < length(labels)
        return
    end
    for cnt = 1:length(labels)
        if angle == 0
            set(text_h(end-cnt+1),'Rotation',0, ...
                'String',labels{cnt},'HorizontalAlignment','center');
        elseif angle > 0
            set(text_h(end-cnt+1),'Rotation',angle, ...
                'String',labels{cnt},'HorizontalAlignment','right');
        elseif angle < 0
            set(text_h(end-cnt+1),'Rotation',angle, ...
                'String',labels{cnt},'HorizontalAlignment','left');
        end
    end
end