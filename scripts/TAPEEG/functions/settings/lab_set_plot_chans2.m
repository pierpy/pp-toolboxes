function [settings,data,skipprocessing] = lab_set_plot_chans(settings,data)

skipprocessing = 0;

if isempty(data)
    skipprocessing = 1;
    return
end

maxdata = size(data,2);
for datanr = 1:maxdata
    PLOT(datanr).Name = ['Data ' num2str(datanr)];
    FORMAT{strcmp(fieldnames(PLOT),'Name')} = 'text';
    PLOT(datanr).Inverse = false;
    if min(data(:,datanr)) < 0
        if min(data(:,datanr)) > -1
            PLOT(datanr).MinValue = -max(abs(data(:,datanr)));
        else
            PLOT(datanr).MinValue = -ceil(max(abs(data(:,datanr))));
        end
    else
        PLOT(datanr).MinValue = min(data(:,datanr));
    end
    if max(abs(data(:,datanr))) < 1
        PLOT(datanr).MaxValue = max(abs(data(:,datanr)));
    else
        PLOT(datanr).MaxValue = ceil(max(abs(data(:,datanr))));
    end
    PLOT(datanr).Threshold1 = 0;
    PLOT(datanr).Threshold2 = 0;
    if min(data(:,datanr)) < 0
        PLOT(datanr).Color1 = 'bluered';
        PLOT(datanr).Color2 = [1 1 1];
        FORMAT{strcmp(fieldnames(PLOT),'Color1')} = {'bluered','autumn','bone','colorcube','cool','copper', ...
            'gray','hot','hsv','jet','iautumn','ibone','icolorcube','icool','icopper','igray','ihot','ihsv','ijet'};
    else
        PLOT(datanr).ColorMode = 'color';
        FORMAT{strcmp(fieldnames(PLOT),'ColorMode')} = {'color','bluered','autumn','bone','colorcube','cool', ...
            'copper','gray','hot','hsv','jet','iautumn','ibone','icolorcube','icool','icopper','igray','ihot','ihsv','ijet'};
        PLOT(datanr).Color1 = [1 0 0];
        PLOT(datanr).Color2 = [1 1 1];
    end
    PLOT(datanr).AddPlot = false;
end

numdiags = ceil(size(PLOT,2) / 20);
if ~exist('FORMAT','var')
    FORMAT = {};
end

for i = 1:numdiags
    Pstart = (i-1)*20+1;
    Pstop = i*20;
    if Pstop > size(PLOT,2)
        Pstop = size(PLOT,2);
    end
    [answer,skipprocessing] = inputsdlg(PLOT(Pstart:Pstop),'Plot settings',FORMAT);
    if skipprocessing == 1
        return
    else
        PLOT(Pstart:Pstop) = answer;
        clearvars answer
    end
end

for datanr = 1:maxdata
    if isfield(PLOT(datanr),'Inverse') & PLOT(datanr).Inverse == true
        if ~isempty(PLOT(datanr).Threshold1)
            tmp = (data(:,datanr) - PLOT(datanr).Threshold1);
            tmp(tmp<0) = 0;
            if ~isempty(PLOT(datanr).Threshold2)
                tmp = tmp / (PLOT(datanr).Threshold2 - PLOT(datanr).Threshold1);
            end
            tmp = tmp/2;
            tmp(tmp>0.5) = 1;
            data(:,datanr) = 1 - tmp;
            clearvars tmp
            PLOT(datanr).MinValue = 0;
            PLOT(datanr).MaxValue = 1;
        elseif max(data(:,datanr)) <= 1 & min(data(:,datanr)) >= 0
            data(:,datanr) = 1 - data(:,datanr);
            PLOT(datanr).MinValue = 0;
            PLOT(datanr).MaxValue = 1;
        elseif min(data(:,datanr)) < 0
            data(:,datanr) = - data(:,datanr);
            tmp = PLOT(datanr).MaxValue;
            PLOT(datanr).MaxValue = -PLOT(datanr).MinValue;
            PLOT(datanr).MinValue = -tmp;
            clearvars tmp
        else
            data(:,datanr) = data(:,datanr).^-1;
            PLOT(datanr).MinValue = 0;
            PLOT(datanr).MaxValue = round(max(data(:,datanr)));
        end
    end
end

plotlist{1} = 1;
for i = 2:maxdata
    if PLOT(i).AddPlot == true
        plotlist{i} = [plotlist{i-1} i];
        plotlist{i-1} = [];
    else
        plotlist{i} = i;
    end
end
tmp = ~cellfun('isempty',plotlist);
plotlist = plotlist(tmp);
clearvars tmp
plotnr2 = plotlist{1};
if length(plotnr2) > 1 & isfield(PLOT(plotnr2(1)),'Threshold1') & (PLOT(plotnr2(1)).Threshold1 > 0 | PLOT(plotnr2(2)).Threshold1 > 0)
    for j = plotnr2
        if isempty(PLOT(j).Threshold2)
            PLOT(j).Threshold2 = 0;
        end
        data(:,j) = (1-data(:,j))*2*abs(PLOT(j).Threshold2-PLOT(j).Threshold1) + PLOT(j).Threshold1;
    end
    settings.pvalues = [PLOT(plotnr2(1)).Threshold2 PLOT(plotnr2(1)).Threshold1];
    settings.pvaluesc = [PLOT(plotnr2(2)).Threshold2 PLOT(plotnr2(2)).Threshold1];
    plotnr2 = plotnr2(1:2);
    if isfield(PLOT,'Color2')
        settings.colors = [{PLOT(plotnr2(1)).Color2} {PLOT(plotnr2(1)).Color1} {PLOT(plotnr2(2)).Color2} {PLOT(plotnr2(2)).Color1}];
    end
elseif isfield(PLOT(plotnr2(1)),'Threshold1') & PLOT(plotnr2(1)).Threshold1 > 0
    plotnr2 = plotnr2(1);
    for j = plotnr2
        if isempty(PLOT(j).Threshold2)
            PLOT(j).Threshold2 = 0;
        end
        data(:,j) = (1-data(:,j))*2*abs(PLOT(j).Threshold2-PLOT(j).Threshold1) + PLOT(j).Threshold1;
    end
    settings.pvalues = [PLOT(plotnr2).Threshold2 PLOT(plotnr2).Threshold1];
    if isfield(PLOT,'Color2')
        settings.colors = [{PLOT(plotnr2(1)).Color2} {PLOT(plotnr2(1)).Color1} {PLOT(plotnr2(1)).Color2} {PLOT(plotnr2(1)).Color1}];
    end
elseif abs(PLOT(plotnr2(1)).MaxValue) > 0
    settings.values = [PLOT(plotnr2(1)).MinValue PLOT(plotnr2(1)).MaxValue];
end

end