% Helper file for lab_plot_IS
%
% Written by F. Hatz 2013

function Chandle = lab_plot_IS_legend(figure1,cfg,mycmap,datanr)

if length(datanr) > 1
    datanr2 = datanr(end);
    datanr = datanr(end-1);
else
    datanr2 = [];
end

if cfg.PLOT(datanr).MaxValue > 0 & cfg.PLOT(datanr).MaxValue > cfg.PLOT(datanr).MinValue
    if isfield(cfg.PLOT,'Threshold1') & cfg.PLOT(datanr).Threshold1 > 0
        if isfield(cfg.PLOT,'Threshold2') & cfg.PLOT(datanr).Threshold2 > 0
            hold on
            circle(0.75,0.044,0.01,mycmap{datanr}(end,:));
            hold on
            circle(0.86,0.044,0.01,mycmap{datanr}(round(size(mycmap{datanr},1)/2),:));
            t1 = text(0.765,0.044,['p <= ' num2str(cfg.PLOT(datanr).Threshold1)],'FontSize',7, ...
                'FontUnit','normalized');
            t2 = text(0.875,0.044,['p <= ' num2str(cfg.PLOT(datanr).Threshold2)],'FontSize',7, ...
                'FontUnit','normalized');
        else
            hold on
            circle(0.86,0.044,0.01,mycmap{datanr}(end,:));
            t1 = text(0.875,0.044,['p <= ' num2str(cfg.PLOT(datanr).Threshold1)],'FontSize',7, ...
                'FontUnit','normalized');
        end
    else
        set(figure1,'Colormap',mycmap{datanr});
        caxis([cfg.PLOT(datanr).MinValue cfg.PLOT(datanr).MaxValue])
        Chandle(1) = colorbar('peer',gca,'location','southoutside','position',[0.732 0.044 0.209 0.05], ...
            'outerposition',[0.697 0.008 0.27 0.06],'FontSize',7);
    end
    if max(cfg.backgroundcolor) == 0
        set(Chandle(1),'XColor',cfg.textcolor);
    end
end
if ~isempty(datanr2)
    if cfg.PLOT(datanr2).MaxValue > 0 & cfg.PLOT(datanr2).MaxValue > cfg.PLOT(datanr2).MinValue
        if isfield(cfg.PLOT,'Threshold1') & cfg.PLOT(datanr2).Threshold1 > 0
            if isfield(cfg.PLOT,'Threshold2') & cfg.PLOT(datanr2).Threshold2 > 0
                hold on
                circle(0.05,0.044,0.01,mycmap{datanr2}(end,:));
                hold on
                circle(0.16,0.044,0.01,mycmap{datanr2}(round(size(mycmap{datanr2},1)/2),:));
                t1 = text(0.065,0.044,['p <= ' num2str(cfg.PLOT(datanr2).Threshold1)],'FontSize',7, ...
                    'FontUnit','normalized');
                t2 = text(0.175,0.044,['p <= ' num2str(cfg.PLOT(datanr2).Threshold2)],'FontSize',7, ...
                    'FontUnit','normalized');
            else
                hold on
                circle(0.05,0.044,0.01,mycmap{datanr2}(end,:));
                t1 = text(0.065,0.044,['p <= ' num2str(cfg.PLOT(datanr2).Threshold1)],'FontSize',7, ...
                    'FontUnit','normalized');
            end
        else
            cbfreeze;
            set(figure1,'Colormap',mycmap{datanr2});
            Chandle(2) = colorbar('peer',gca,'location','southoutside','position',[0.045 0.044 0.209 0.05], ...
                'outerposition',[0.045 0.008 0.27 0.06],'FontSize',7);
        end
        if max(cfg.backgroundcolor) == 0
            set(Chandle(2),'XColor',cfg.textcolor);
        end
    end
end
if exist('Chandle','var')
    U = get(gcf,'UserData');
    U.Hcolorbar = Chandle;
    set(gcf,'UserData',U);
end

function circle(x,y,r,ccolor)

hold on
th = 0:pi/180:2*pi;
xunit = r * cos(th) + x;
yunit = r * sin(th) + y;
fill(xunit, yunit,ccolor,'EdgeColor',ccolor);
hold off