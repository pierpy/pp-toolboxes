function RaguShowMSFit(h)

% Copyright 2009-2011 Thomas Koenig
% distributed under the terms of the GNU AFFERO General Public License

out = get(h,'UserData');
figure(h);

if (isfield(out,'strF1'))
    strF1 = out.strF1;
else
    strF1 = 'Factor 1';
end

if (isfield(out,'strF2'))
    strF2 = out.strF2;
else
    strF2 = 'Factor 2';
end

handles = guidata(gcf);
handles.CurrentView = 'MSFit';
guidata(gcf,handles);


factors = {'Main effect' strF1 strF2 [strF1 ' * ' strF2]};

cm = MyColors;

nMaps = size(out.MSMaps,1);
dx = 0.7 / nMaps;
for i = 1:size(out.MSMaps,1)
    h = subplot('Position',[(i-1)*dx + 0.31,0.85,dx-0.02,0.15]);
    RaguDSPMap(out.MSMaps(i,:),out.Channel,out.MapStyle,'NoScale','Resolution',3,'Title',sprintf('Class %i',i));
%    set(get(h,'Children'),'ButtonDownFcn',{@MapBtCallback,out.MSMaps(i,:),out.Channel,out.MapStyle,cm(i,:),i});
%    title(sprintf('Class %i',i),'Units','normalized');
    axis off
    h2 = uicontrol('Units','normalized','Position',[(i-1)*dx + 0.31,0.81,dx-0.02,0.035],'String',sprintf('Select class %i',i),'BackgroundColor',cm(i,:),'Tag','MSStatsTable');
%    subplot('Position',[(i-1)*dx + 0.31,0.82,dx-0.02,0.018]);
%    h2 = patch([0 1 1 0 0],[0 0 1 1 0],cm(i,:));
%    axis([0 1 0 1]);
%    set(h2,'Linewidth',1);
%    axis off
    if (isfield(out,'pMSStats'))
        if ~isempty(out.pMSStats)
            ud.Stats = squeeze(out.pMSStats(:,:,i,1,:));
        else
            ud.Stats = [];
        end
    else
        ud.Stats = [];
    end
    ud.Class = i;
    ud.Factors = factors;
    set(h2,'UserData',ud,'Callback',{@ShowStats,gcf,i});
end
ShowAssignment();


function MapBtCallback(obj, event, data,montage,style,col,index)

figure;
subplot('Position',[0.1 0.2 0.8 0.8]);

RaguDSPMap(data,montage,style,'NoScale','Resolution',1);%,'Label',lbl');
subplot('Position',[0.1 0.1 0.8 0.05]);
title(sprintf('Class %i',index),'Units','normalized');
h2 = patch([0 1 1 0 0],[0 0 1 1 0],col);
axis([0 1 0 1]);
set(h2,'Linewidth',1);
axis off

function cm = MyColors()
    cm = [...
0.00 0.00 1.00;...
0.00 0.50 0.00;...
1.00 0.00 0.00;...
0.00 0.75 0.75;...
0.75 0.00 0.75;...
0.75 0.75 0.00;...
0.25 0.25 0.25;...
0.00 0.00 0.50;...
1.00 1.00 0.00;...
0.50 0.50 1.00;...
1.00 0.50 0.50;...
0.50 1.00 0.50;...
0.75 0.50 0.25;...
0.00 1.00 1.00;...
1.00 0.00 1.00;...
0.75 0.75 0.25;...
1.00 0.75 0.25;...
0.50 0.00 0.00;...
0.50 0.50 0.50;...
0.00 0.50 0.50...
            ];


function ShowStats(hObject,eventdata,handles,MSClass)
    
ud = get(hObject,'UserData');

if ~isempty(ud.Stats)
    onset    = squeeze(ud.Stats(:,:,1))';
    offset   = squeeze(ud.Stats(:,:,2))';
    duration = squeeze(ud.Stats(:,:,3))';
    auc      = squeeze(ud.Stats(:,:,4))';
    cog      = squeeze(ud.Stats(:,:,5))';
    gfp      = squeeze(ud.Stats(:,:,6))';
else
    onset    = repmat({'Show'},4,2);
    offset   = repmat({'Show'},4,2);
    duration = repmat({'Show'},4,2);
    auc      = repmat({'Show'},4,2);
    cog      =repmat({'Show'},4,2);
    gfp      = repmat({'Show'},4,2);
end

PlaceMyTable('Onset'            ,1,onset   ,ud);
PlaceMyTable('Offset'           ,2,offset  ,ud);
PlaceMyTable('Duration'         ,3,duration,ud);
PlaceMyTable('Area under curve' ,4,auc     ,ud);
PlaceMyTable('Center of gravity',5,cog     ,ud);
PlaceMyTable('Mean GFP'         ,6,gfp     ,ud);
ShowAssignment(MSClass);


function PlaceMyTable(txt,pos,data,ud)

h = findobj('tag','MSStatsTable','UserData',pos);
for i = 1:numel(h)
    delete(h(i));
end
mydata = cell(4,3);
for i = 1:4
    mydata(i,1) = ud.Factors(i);
    if isnumeric(data(i,1))
        mydata(i,2) = {data(i,1)};
        mydata(i,3) = {data(i,2)};
    else
        mydata(i,2) = data(i,1);
        mydata(i,3) = data(i,2);
    end
end

dx = 0.96 / 6;
sx = 0.96 / 6.1;

%uitable('Units','pixel','Position',[30 (6-pos) *110 + 30 265 90],'Data',mydata,'ColumnFormat',{'char' 'numeric' 'numeric'},'ColumnName',{txt,'Main','Group'},'RowName',[],'tag','MSStatsTable','CellSelectionCallback',{@CellClicked,gcf,ud,pos},'TooltipString','Click on cell to display','Visible','on','ColumnWidth',{120 70 70},'UserData',pos);
uitable('Units',' normalized','Position',[0.01 0.02 + (6-pos) * dx 0.27 sx],'Data',mydata,'ColumnFormat',{'char' 'numeric' 'numeric'},'ColumnName',{txt,'Main','Group'},'RowName',[],'tag','MSStatsTable','CellSelectionCallback',{@CellClicked,gcf,ud,pos},'TooltipString','Click on cell to display','Visible','on','ColumnWidth',{120 70 70},'UserData',pos);
function UnselectTable(pos)

h = findobj('tag','MSStatsTable','UserData',pos);

for i = 1:numel(h)
    data = get(h,'Data');
    cb = get(h,'CellSelectionCallback');
    set(h,'Data',[],'CellSelectionCallback',[]);
    set(h,'Data',data,'CellSelectionCallback',cb);
    
end




function CellClicked(hObject,eventdata,handles,ud,tbl)

if ~isempty(ud.Stats)
    onset    = squeeze(ud.Stats(:,:,1))';
    offset   = squeeze(ud.Stats(:,:,2))';
    duration = squeeze(ud.Stats(:,:,3))';
    auc      = squeeze(ud.Stats(:,:,4))';
    cog      = squeeze(ud.Stats(:,:,5))';
    gfp      = squeeze(ud.Stats(:,:,6))';
else
    onset    = repmat({'Show'},4,2);
    offset   = repmat({'Show'},4,2);
    duration = repmat({'Show'},4,2);
    auc      = repmat({'Show'},4,2);
    cog      =repmat({'Show'},4,2);
    gfp      = repmat({'Show'},4,2);
end


for i = 1:6
    if (tbl ~= i)
        UnselectTable(i);
    end
end
% if (tbl ~= 1)
%     PlaceMyTable('Onset',1,onset,ud);
%  end
%  
%  if (tbl ~= 2)
%     PlaceMyTable('Offset',2,offset,ud);
% end
% 
% if (tbl ~= 3)
%     PlaceMyTable('Duration',3,duration,ud);
% end
% 
% if (tbl ~= 4)
%     PlaceMyTable('Area under curve' ,4,auc,ud);
% end
% 
% if (tbl ~= 5)
%     PlaceMyTable('Center of gravity',5,cog,ud);
% end
% 
% if (tbl ~= 6)
%     PlaceMyTable('Mean GFP',6,gfp,ud);
% end

ShowAssignment(ud.Class,eventdata.Indices,tbl)


function ShowAssignment(SelectedMSClass,Indices,tbl)

if nargin > 1
    if(size(Indices,1) == 0)
        return;
    end
end

out = get(gcf,'UserData');
Group = out.IndFeature;
GroupLabels = out.GroupLabels;

SelDesign = out.Design;
SelDesign(isnan(SelDesign(:,1)),:) = [];

DoF1    = (numel(unique(SelDesign(:,1)))> 1);
DoF2 = out.TwoFactors;

if nargin < 2
    if (numel(unique(Group(~isnan(Group)))) > 1)
        Indices(2) = 3;
    else
        Indices(2) = 1;
    end
   
    Indices(1) = DoF1+2*DoF2+1;
    
end

if Indices(2) < 3 
    Group(~isnan(Group)) = 1;
end

if nargin > 1 && Indices(2) < 3
    GroupLabels = [];
end

AllC = ones(size(out.Design,1),1);
AllC(isnan(out.Design(:,1)),1) = nan;

switch(Indices(1))
    case 1
        [gm,outlbl] = RaguGrandMeans(out.V,Group,AllC,out.DLabels1,out.DLabels2,GroupLabels ); % TK 14.3.2013
    case 2
        [gm,outlbl] = RaguGrandMeans(out.V,Group,out.Design(:,1),out.DLabels1,out.DLabels2,GroupLabels );

    case 3
        [gm,outlbl] = RaguGrandMeans(out.V,Group,out.Design(:,2),out.DLabels2,out.DLabels1,GroupLabels );
    case 4
        [gm,outlbl] = RaguGrandMeans(out.V,Group,out.Design,out.DLabels1,out.DLabels2,GroupLabels );
end

nc = size(gm,2);
ng = size(gm,1);

if out.DoAnTopFFT == 0
    lbltime = (0:(size(gm,4)-1)) * out.DeltaX + out.TimeOnset;
else
    lbltime = (0:(size(gm,4)-1)) * out.DeltaX;
end

xlim = [min(lbltime) max(lbltime)];
gfp = std(gm,1,3);

mx = max(gfp(:));

cm = vga;
cm = cm(2:end,:);

[MSClass,MSFit] = RaguFitMicrostates(gm,out.MSMaps,out.bSmoothLabels,out.nWindowSize,out.LabelPenalty,out.StartFrame,out.EndFrame);
[on,off,dur,auc,cog,msgfp] = RaguMSOnOffsets(MSClass,size(out.MSMaps,1),MSFit,gfp, 0);

infotxt = '';

dx =  0.68 / ng;
dy = 0.77 / nc;

for c = 1:nc
    for g = 1:ng

%        sph = subplot('Position',[(0.32 + (g-1) * dx) (0.85 - (c) * dy)  (dx - 0.04) (dy-0.04)]);
        sph = subplot('Position',[(0.32 + (g-1) * dx) (0.81 - (c) * dy)  (dx - 0.02) (dy-0.03)]);
        if c == 1 && g == 1
           plot(1,mx);
           ylim = get(gca,'YLim');
           ylim(1) = 0;
           cla
        end
    
        d = zeros(size(out.MSMaps,1),size(MSClass,3));

        for t = 1:size(d,2)
            if MSClass(g,c,t) > 0
                d(MSClass(g,c,t),t) = MSFit(g,c,t);
            end
        end

        plot(lbltime,squeeze(gfp(g,c,:)),'-k','Linewidth',1);
        hold on
                
            
        bh = bar(lbltime,d',1,'stacked');
                
        bcm = MyColors;
    
        for i = 1:numel(bh)
            set(bh(i),'FaceColor',bcm(i,:),'LineStyle','none','ButtonDownFcn',{@GraphBtCallback,lbltime,d,bcm,outlbl{g,c},xlim,ylim,size(cm,1)});
        end
        
        set(sph,'ButtonDownFcn',{@GraphBtCallback,lbltime,d,bcm,outlbl{g,c},xlim,ylim,size(cm,1)});
        
%        shading flat
%        colormap(cm);
        set(gca,'YLim',ylim,'XLim',xlim,'CLim',[1 size(cm,1)],'TickDir','in','Box','off');
        xlabel(out.txtX);
        ylabel('GFP');
        
%        if(g > 1)
%            set(gca,'YTickLabel',[]);
%        end

%        if(c < nc)
%            set(gca,'XTickLabel',[]);
%        end
       
        if nargin > 2
            switch tbl
                case 1
                    infotxt = sprintf(' (Onset = %1.1f %s)',on(g,c,SelectedMSClass)* out.DeltaX + out.TimeOnset,out.txtX);
                case 2
                    infotxt = sprintf(' (Offset = %1.1f  %s)',off(g,c,SelectedMSClass)* out.DeltaX + out.TimeOnset,out.txtX);
                case 3
                    infotxt = sprintf(' (Duration = %1.1f  %s)',dur(g,c,SelectedMSClass)* out.DeltaX,out.txtX);
                case 4
                    infotxt = sprintf(' (AUC = %1.1f  %s x uV)',auc(g,c,SelectedMSClass)* out.DeltaX,out.txtX);
                case 5
                    infotxt = sprintf(' (COG = %1.1f  %s)',cog(g,c,SelectedMSClass)* out.DeltaX + out.TimeOnset,out.txtX);
                case 6
                    infotxt = sprintf(' (GFP = %1.4f  uV)',msgfp(g,c,SelectedMSClass));

            end
        end
        text(mean(xlim),ylim(2),[outlbl{g,c} infotxt],'Interpreter','none','VerticalAlignment','top','HorizontalAlignment','center');
%        title([outlbl{g,c} infotxt],'Interpreter','none');
                
        if nargin > 0
            dsptime = lbltime;
            for t = 1:(numel(dsptime)-1)
                if (d(SelectedMSClass,t) == 0) && (d(SelectedMSClass,t+1) > 0)
                    dsptime(t+1) = (dsptime(t+1) + dsptime(t))/2;
                    dsptime(t  ) = dsptime(t+1);
                end
                if (d(SelectedMSClass,t) > 0) && (d(SelectedMSClass,t+1) == 0)
                    dsptime(t  ) = (dsptime(t+1) + dsptime(t))/2;
                    dsptime(t+1) = dsptime(t);
                end
            end
                
            plot(dsptime,d(SelectedMSClass,:),'-k','Linewidth',3);

        end
                
%        xl = get(gca,'XTick');
%        dlt = xl(2) - xl(1);
%        xnl = (xl+1):xl:numel(time);
%        set(gca,'XTickLabel',num2str(lbltime(xnl)'),'XTick',xnl);
        hold off
        freezeColors();
    end
end

function GraphBtCallback(obj, event, lbltime,d,bcm,tit,xlim,ylim,clim)

figure

bh = bar(lbltime,d',1,'stacked');
                
for i = 1:numel(bh)
    set(bh(i),'FaceColor',bcm(i,:),'LineStyle','none');
end

title(tit);

set(gca,'YLim',ylim,'XLim',xlim,'CLim',[1 clim]);