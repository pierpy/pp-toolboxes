function ShowAnTopographyResults(out,fig)

% Copyright 2009-2011 Thomas Koenig
% distributed under the terms of the GNU AFFERO General Public License

handles = guidata(gcf);
handles.CurrentView = 'TCT';
guidata(gcf,handles);

if out.MeanInterval == 0
    time = ((out.StartFrame:out.EndFrame)-1) * out.DeltaX + out.TimeOnset;
else
    time = 1;
end
mx = max(out.MeanGFP(:));

mx = mx * 1.2;
bc = 0.2;
bcol = [bc bc bc];

if nargin < 2
    figure(101);
    clf;
else
    figure(fig);
end

if ~out.ContBetween
    b = unique(out.IndFeature(~isnan(out.IndFeature)));
    for i = 1:numel(b)
        Group{i} = find(out.IndFeature == b(i));
    end
else
    Group{1} = 1:numel(out.IndFeature(~isnan(out.IndFeature)));
end 

VMean = nan(numel(Group),size(out.V,2),size(out.V,3),numel(time));
VTVal = nan(numel(Group),size(out.V,2),size(out.V,3),numel(time));

for c = 1:size(out.V,2)
    for grp = 1:numel(Group)
        if out.MeanInterval == 0
            VMean(grp,c,:,:) = mean(out.V(Group{grp},c,:,out.StartFrame:out.EndFrame),1);
            VTVal(grp,c,:,:) = tvalue(out.V(Group{grp},c,:,out.StartFrame:out.EndFrame),[],1);
        else
            VMean(grp,c,:,1) = mean(mean(out.V(Group{grp},c,:,out.StartFrame:out.EndFrame),4),1);
            VTVal(grp,c,:,1) = tvalue(mean(out.V(Group{grp},c,:,out.StartFrame:out.EndFrame),4),[],1);
        end
    end
end



nGroups = size(out.MeanGFP,1);
nConds  = size(out.MeanGFP,2);

dx = 1/nGroups;
dy = 1/nConds;

if (out.MeanInterval == 0)
    ButtonName = questdlg('Apply duration threshold?','Correction for multiple testing','Yes', 'No', 'No');
    switch ButtonName
        case 'Yes'
            CritDur = out.TCTCritDuration;
        case 'No'
            CritDur = ones(size(out.TCTCritDuration)) * 10000;
    end
else
    CritDur = ones(size(out.TCTCritDuration)) * 10000;
end

rsig = zeros(nGroups,nConds,size(out.MeanGFP,3));

if out.DoAnTopFFT == 0
    for g = 1:nGroups
        for c = 1:nConds
            sig = squeeze((out.pTopCons(g,c,:)) > out.Threshold);
            rsig(g,c,:) = zeros(size(sig));
            isOn = false;
            for tm = 1:numel(sig);
                if sig(tm) == false
                    if isOn == false
                        tStart = tm;
                    end
                        isOn = true;
                else
                    if isOn == true && (tm - tStart - 1) >= CritDur(g,c)
                        rsig(g,c,tStart:(tm-1)) = 1;
                    end
                    isOn = false;
                end
            end
            if isOn == true && (tm - tStart) >= CritDur(g,c)
                rsig(g,c,tStart:tm) = 1;
            end
        end
    end
end

for c = 1:nConds
    for g = 1:nGroups
        xp = (g-1) * dx;
        yp = 1-c*dy;
        ah2 = subplot('Position',[xp+0.75*dx,yp+0.2*dy,dx*0.2,dy * 0.6]);
        ah  = subplot('Position',[xp+0.1*dx,yp+0.2*dy,dx*0.6,dy * 0.6]);
        axis(ah2,'off');
        
        set(ah,'Tag','AnTopTextBox');
        IsSig = squeeze(rsig(g,c,:));
        pValue = squeeze(out.pTopCons(g,c,:));
        IsSig = IsSig - IsSig .* pValue;
        ToPlot = [pValue,IsSig];
        Rest = 1-ToPlot;
        ToPlot = [ToPlot Rest];
        if out.MeanInterval == 0
            h = bar(time',ToPlot * mx,'stacked');
        else
            h = bar([ToPlot;nan(size(ToPlot))] * mx,'stacked');
            xlim([0.5 1.5]);
        end
        
        PlotInfo.EEG = squeeze(VMean(g,c,:,:));
        PlotInfo.time = time;
        PlotInfo.Channel = out.Channel;
        PlotInfo.mx = mx;
        PlotInfo.ax = ah2;
        PlotInfo.MapStyle = out.MapStyle;
        PlotInfo.Label = out.txtX;
        PlotInfo.p = pValue;
        PlotInfo.tValues = squeeze(VTVal(g,c,:,:));
%        set(ah,'ButtonDownFcn',{@ShowAnTopMap,ah,PlotInfo});
                
        gray = [0.7 0.7 0.7];
        set(h(1),'BarWidth',1,'EdgeColor',gray,'FaceColor',gray);
        set(h(2),'BarWidth',1,'EdgeColor',[0 0 0],'FaceColor',[0 1 0],'LineStyle','none');
        set(h(3),'BarWidth',1,'EdgeColor',[0 0 0],'FaceColor',[1 1 1],'LineStyle','none');

        hold on
        plot([(time(1)-0.5) (max(time)+0.5)],[out.Threshold out.Threshold] * mx,'-r');
        plot(time,squeeze(out.MeanGFP(g,c,:)),'LineWidth',1,'Color',[0,0,0]);
        axis([(time(1)-0.5) (max(time)+0.5) 0 mx]);
        
        if out.ContBetween == false
            s = [char(out.GroupLabels{b(g),1}) ': ' char(out.conds{c,1}) '  Overall p: ' num2str(out.TCTPHitCount(g,c))];
        else
            s = [char(out.conds{c,1}) '  Overall p: ' num2str(out.TCTPHitCount(g,c))];
        end
        if out.DoAnTopFFT == 0
            h = title([s ' (min duration: ' num2str(out.TCTCritDuration(g,c) * out.DeltaX) ' ' out.txtX ')']);
        else
            h = title(s);
        end
        set(h,'Interpreter','none','Fontsize',8);
        xlabel(out.txtX,'Tag','AnTopTextBox');
        ylabel('GFP','Tag','AnTopTextBox');
        if out.MeanInterval == 0
            ahc = get(gca,'Children');
            for kid = 1:numel(ahc)
                set(ahc(kid),'ButtonDownFcn',{@ShowAnTopMap,ah,PlotInfo});
            end
        else        
            ShowAnTopMapIdx(PlotInfo,1);
        end
    end
end

if nConds > 1
    for c = 1:(nConds-1)
        yp = 1-c*dy;
        annotation('line',[0 1],[yp yp],'Tag','AnTopTextBox');
    end
end

if nGroups > 1
    for g = 1:(nGroups-1)
        xp = g * dx;
        annotation('line',[xp xp],[0 1],'Tag','AnTopTextBox');
    end
end
   

function ShowAnTopMap(hObject, eventdata, ax,Info)

ch = get(ax,'UserData');

if ishandle(ch)
    delete(ch);
end

p = get(ax,'CurrentPoint');

[~,idx] = min(abs(Info.time - p(1,1)));

y = get(ax,'YLim');

ch = line([p(1,1) p(1,1)],y,'Color',[0 0 1],'LineWidth',1);
set(ax,'UserData',ch);

ShowAnTopMapIdx(Info,idx);

function ShowAnTopMapIdx(Info,idx)

h = gca;
axes(Info.ax);

RaguDSPMap(Info.EEG(:,idx),Info.Channel,Info.MapStyle,'Step',Info.mx/4,'Title',sprintf('TCT: %3.2f %s, p = %4.4f',Info.time(idx),Info.Label,Info.p(idx)),'tValue',Info.tValues(:,idx));
title(sprintf('%3.2f %s, p = %4.4f',Info.time(idx),Info.Label,Info.p(idx)),'Interpreter','none','Fontsize',8);
axes(h);




