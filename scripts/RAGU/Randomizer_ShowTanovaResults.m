function Randomizer_ShowTanovaResults(out,h)

% Copyright 2009-2011 Thomas Koenig
% distributed under the terms of the GNU AFFERO General Public License

handles = guidata(gcf);
handles.CurrentView = 'TANOVA';
guidata(gcf,handles);

if out.DoGFP == 0
    pTanova = out.PTanova;
else
    pTanova = out.GFPPTanova;
end

time = (0:(size(out.V,4)-1)) * out.DeltaX + out.TimeOnset;
axislabel = out.txtX;

axismin = time(out.StartFrame);
axismax = time(out.EndFrame);

Threshold = ones(2,8) * out.Threshold;

if (isfield(out,'CritFDR_p') && out.DoGFP == 0)
    if strcmp(questdlg('Use FDR Threshold?','show TANOVA results','Yes','No','Yes'),'Yes')
        Threshold = out.CritFDR_p;
    end
end    

if nargin < 2
    h = figure(103);
    set(h,'Units','normalized','Position',[0.05 0.05 0.9 0.8],'MenuBar','figure','ToolBar','figure');
    clf

    mymenu = uimenu('Label','Export');
    uimenu(mymenu,'Label','Save as Metafile','Callback',{@SaveFigure,'-dmeta','*.wmf', 'Save figure to Metafile'});
    uimenu(mymenu,'Label','Save as Bitmap','Callback',{@SaveFigure,'-dbitmap','*.bmp', 'Save figure to bitmap'});
    uimenu(mymenu,'Label','Copy as Metafile','Callback',{@SaveFigure,'-dmeta',[], []});
    uimenu(mymenu,'Label','Copy as Bitmap','Callback',{@SaveFigure,'-dbitmap',[], []});
    uimenu(mymenu,'Label','Export to text','Callback',{@ExportFigure,'*.txt', 'Save as textfile',PTanova});
end

bc = 0.7;
bcol = [bc bc bc];

SelDesign = out.Design;
SelDesign(isnan(SelDesign(:,1)),:) = [];

DoF1    = (numel(unique(SelDesign(:,1)))> 1);
DoF2    = (numel(unique(SelDesign(:,2)))> 1);
if out.TwoFactors == 0
    DoF2 = 0;
end

if DoF1 && DoF2
    nc = 4;
elseif DoF1 && ~DoF2
    nc = 2;
else
    nc = 1;
end

if numel(unique(out.IndFeature(~isnan(out.IndFeature))))> 1
    ng = 2;
else
    ng = 1;
end

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

t = {'' [out.IndName ' main effect'] strF1 [strF1 ' * ' out.IndName] strF2 [strF2 ' * ' out.IndName] [strF1 ' * ' strF2] [strF1 ' * ' strF2 ' * ' out.IndName]};

hp = zeros(2,4);

d = get(gcf,'UserData');
if ng == 2
    d.MainGroupLabel = 'All';
else
    b = unique(out.IndFeature(~isnan(out.IndFeature)));
    d.MainGroupLabel = char(out.GroupLabels{b(1)});
end

DataToUse = find(~isnan(out.Design(:,1)));
gfp = sqrt(mean(out.V(:,DataToUse,:,:).^2,3));

% cmeans: Covariance maps for each group
% gmeans: Within effects across groups
% means: Means within group and condition

if out.ContBetween == false
    b = unique(out.IndFeature(~isnan(out.IndFeature)));
    d.GLabels = out.GroupLabels(b);
    nGroups = sum(~isnan(b));
    
    d.means    = zeros(nGroups,numel(DataToUse),size(out.V,3),size(out.V,4));
    d.meansGFP = zeros(nGroups,numel(DataToUse),1,size(out.V,4));
    
    if d.ContF1 
        d.cmeans    = zeros(nGroups,2,size(out.V,3),size(out.V,4));
        d.cmeansGFP = zeros(nGroups,1,1            ,size(out.V,4));
    end
    
    % These are the within effects for each group
    for g = 1:nGroups
        if d.ContF1
            [cm,cmgfp] = MakeCovMapForDisplay(out.V(out.IndFeature == b(g),DataToUse,:,:),out.Design(DataToUse,1),[]);
            
            d.cmeans(g,:,:,:)   = cm;
            d.cmeansGFP(g,:,:,:) = cmgfp;
        end
        d.means(g,:,:,:)    = mean(out.V(out.IndFeature == b(g),DataToUse,:,:),1);
        d.meansGFP(g,:,:,:) = mean(  gfp(out.IndFeature == b(g),:        ,:,:),1);
    end
else
    % These are the within effects for each group
    if d.ContF1
        d.cmeans    = zeros(2,2,size(out.V,3),size(out.V,4));
        d.cmeansGFP = zeros(2,2,1,            size(out.V,4));
        
        [cm,cmgfp] = MakeCovMapForDisplay(out.V(~isnan(out.IndFeature),DataToUse,:,:),out.Design(DataToUse,1),out.IndFeature(~isnan(out.IndFeature)));
        
        d.cmeans    = cm;
        d.cmeansGFP = cmgfp;
    end
%    d.means    = zeros(1,numel(DataToUse),size(out.V,3),size(out.V,4));
%    d.meansGFP = zeros(1,numel(DataToUse),1            ,size(out.V,4));
    d.means    = zeros(2,numel(DataToUse),size(out.V,3),size(out.V,4));
    d.meansGFP = zeros(2,numel(DataToUse),1            ,size(out.V,4));

    for c = 1:numel(DataToUse)
        [cm,cmgfp] = MakeCovMapForDisplay(out.V(~isnan(out.IndFeature),DataToUse(c),:,:),[],out.IndFeature(~isnan(out.IndFeature)));
        d.means(:,c,:,:) = cm;
        d.meansGFP(:,c,:,:) = cmgfp;
    end
end

% These are the within effects across groups

if d.ContF1
    [cm,cmgfp] = MakeCovMapForDisplay(out.V(~isnan(out.IndFeature),DataToUse,:,:),out.Design(DataToUse,1),[]);
    d.gmeans = cm;
    d.gmeansGFP = cmgfp;
    
else
    d.gmeans    = mean(out.V(~isnan(out.IndFeature),DataToUse,:,:),1);
    d.gmeansGFP = mean(  gfp(~isnan(out.IndFeature),        :,:,:),1);
end
    
    
d.DSPconds = out.conds;
d.DSPconds(isnan(out.Design(:,1))) = [];
d.ng = ng;
d.nc = nc;
d.titles = t;
d.Averaged = 0;
d.SelDesign = SelDesign;
d.ContBetween = out.ContBetween;
d.axislabel = axislabel;
d.time = time;

if size(d.V,4) == 1
    out.MeanInterval = 1;
end

if out.MeanInterval
    d.Averaged = 1;
    
    d.means  = mean( d.means(:,:,:,out.StartFrame:out.EndFrame),4);
    d.gmeans = mean(d.gmeans(:,:,:,out.StartFrame:out.EndFrame),4);
    d.meansGFP = mean(d.meansGFP(:,:,:,out.StartFrame:out.EndFrame),4);
    d.gmeansGFP = mean(d.gmeansGFP(:,:,:,out.StartFrame:out.EndFrame),4);
    
    for i = 1:ng
        for j = 1:nc
            if (i * j > 1) || (nc *ng == 1)
                hp(i,j) = RaguTanovaSubplot(nc,ng,j,i);
                cla(hp(i,j));
                h = bar(0.5,squeeze(pTanova(i,j,1,1)));
                set(h,'BarWidth',1,'EdgeColor',[0 0 0],'FaceColor',[0 0 0]);
                set(h,'Tag','TanovaP');
                title([t{2*j+i-2} ': p= ' sprintf('%3.3f',pTanova(i,j))],'Interpreter','none','FontSize',8);
                axis([0 1 0 1]);
                set(hp(i,j),'Tag','subplot','Userdata',t{2*j+i-2});
                xlabel(d.txtX);
                ylabel('p-Value');
            end
        end
    end
else    
    rsig = zeros(ng,nc,size(pTanova,3));
    for i = 1:ng
        for j = 1:nc
            CritDuration = inf;
            if out.DoGFP == 0
                if isfield(out,'CritDuration')
                    CritDuration = out.CritDuration(i,j);
                end
            else
                if(isfield(out,'CritDurationGFP'))
                    CritDuration = out.CritDurationGFP(i,j);
                end
            end
            sig = squeeze(pTanova(i,j,:,1) > out.Threshold);
            rsig(i,j,:) = zeros(size(sig));
            isOn = false;
            for tm = 1:numel(sig);
                if sig(tm) == false
                    if isOn == false
                        tStart = tm;
                    end
                    isOn = true;
                else
                    if isOn == true && (tm - tStart - 1) >= CritDuration
                        rsig(i,j,tStart:(tm-1)) = 1;
                    end
                    isOn = false;
                end
            end
            if isOn == true && (tm - tStart) >= CritDuration
                rsig(i,j,tStart:tm) = 1;
            end
        end
    end
    for i = 1:ng
        for j = 1:nc
            if (i * j > 1) || (nc *ng == 1)
                hp(i,j) = RaguTanovaSubplot(nc,ng,j,i);
                set(hp(i,j),'Tag','subplot');
                h = bar(time,double(squeeze(pTanova(i,j,:,1) > Threshold(i,j))));
                set(h,'BarWidth',1,'EdgeColor',bcol,'FaceColor',bcol);
                hold on
                h = bar(time,squeeze(rsig(i,j,:)));
                set(h,'BarWidth',1,'EdgeColor',[0 1 0],'FaceColor',[0 1 0]);
                tanh = plot(time,squeeze(pTanova(i,j,:,1)),'-k');
                set(tanh,'Tag','TanovaP');
                title(t{2*j+i-2},'Interpreter','none','FontSize',8);
                axis([axismin axismax 0 1]);
                tc = (axismin + axismax) / 2;
                lh = plot([tc tc], [0 1],'-r');
                set(lh,'Tag','Cursor');
                set(hp(i,j),'Tag','subplot','Userdata',t{2*j+i-2});
                xlabel(d.txtX);
                ylabel('p-Value');

            end
        end
    end
end

if (out.Normalize == 2) && (out.DoGFP == 0)
    d.means = NormDimL2(d.means,3);
    d.gmeans = NormDimL2(d.gmeans,3);
end

mxg =  0;
mxgr = 0;

for i = 1:size(d.means,4)
    maps = reshape(d.means(:,:,:,i),size(d.means,1)*size(d.means,2),size(d.means,3));
    maps = maps - repmat(mean(maps),size(maps,1),1);
    
    maplen = sqrt(sum(maps.*maps,2));
    mxg = max([mxg;maplen]);
end

if(out.ContBetween == true)
    for i = 1:size(d.gmeans,4)
        maps = squeeze(d.gmeans(1,:,:,i));
        maps = maps - repmat(mean(maps),size(maps,1),1);
        maplen = sqrt(sum(maps.*maps,2));
        mxgr = max([mxgr;maplen]);
    end
end
    
    
d.mx = mxg;
d.mxgr = mxgr;
d.mdsScaleFactor = 1;

if isfield(out,'Channel')
    d.Channel = out.Channel;
end

set(gcf,'UserData',d);

if (out.DoGFP == 0)
    uicontrol('Style','pushbutton','String','Zoom','Units','normalized','Position',[0.86 0.21 0.03 0.03],'Callback',{@RescaleMDS,1,true},'Tag','Transient');
    uicontrol('Style','pushbutton','String','+','Units','normalized','Position',[0.86 0.18 0.03 0.03],'Callback',{@RescaleMDS,.75,false},'Tag','Transient');
    uicontrol('Style','pushbutton','String','-','Units','normalized','Position',[0.86 0.15 0.03 0.03],'Callback',{@RescaleMDS,1.5,false},'Tag','Transient');
end

%if (~out.MeanInterval)
    for i = 1:ng
        for j = 1:nc
            if ~strcmp(get(hp(i,j),'Type'),'axes')
                continue
            end
            set(hp(i,j),'ButtonDownFcn',{@ShowTanovaBtCallback, i,j});
            kids = get(hp(i,j),'Children');
            for c = 1:numel(kids)
                set(kids(c),'ButtonDownFcn',{@ShowTanovaBtCallbackKid, i,j});
            end
        end
    end
    
bDone = false;
for i = 1:2
    for j = 1:4
        if(hp(i,j) ~= 0) && bDone == false
            ShowTanovaBtCallback(hp(i,j), 0, i,j);
            bDone = true;
        end
    end
end


function [cm,gfp] = MakeCovMapForDisplay(data,WithinPred,BetweenPred)

    gfp = sqrt(mean(data.^2,3));
    data = data - repmat(mean(mean(data,2),1),[size(data,1) size(data,2) 1 1]); % We center across the conditions and within group
    
    DubWithin = false;
    if ~isempty(WithinPred)
        DubWithin = true;
        WithinPred = WithinPred - mean(WithinPred);
        WithinPred = repmat(WithinPred',[size(data,1),1,1,size(data,4)]);
    
    else
        WithinPred = ones(size(gfp));
    end
    
    DubBetween = false;
    if ~isempty(BetweenPred)
        DubBetween = true;
        BetweenPred = BetweenPred - mean(BetweenPred);
        BetweenPred = repmat(BetweenPred,[1,size(data,2),1,size(data,4)]);
    else
        BetweenPred = ones(size(gfp));
    end
    gfp = mean(mean(gfp.*WithinPred.*BetweenPred,2),1);
    
    WithinPred  = repmat(WithinPred ,[1,1,size(data,3),1]);
    BetweenPred = repmat(BetweenPred,[1,1,size(data,3),1]);
    cm = mean(mean(data.*WithinPred.*BetweenPred,2),1);
    
    if DubBetween == true;
        cm(2,:,:,:) = -cm(1,:,:,:);
        gfp(2,:,:,:) = - gfp(1,:,:,:);
    end
    
    if DubWithin == true;
        cm(:,2,:,:) = -cm(:,1,:,:);
        gfp(:,2,:,:) = -gfp(:,1,:,:);
    end
    
            

function ExportFigure(obj, event,mask,comment,PTanova)

[filename, pathname] = uiputfile(mask,comment);
if isequal(filename,0) || isequal(pathname,0)
    return
end

[fp,err] = fopen(fullfile(pathname, filename),'wt');

if fp == -1
    errordlg(err);
    return
end
    
d = get(gcf,'UserData');

fprintf(fp,'Time');

for i = 1:d.ng
    for j = 1:d.nc
        if (i * j > 1)
            fprintf(fp,'\t%s',d.titles{d.ng*j+i-2});
        end
    end
end


for t = 1:size(PTanova,3);
    fprintf(fp,'\n%5.2f',d.time(t));
    for i = 1:d.ng
        for j = 1:d.nc
            if (i * j > 1)
                fprintf(fp,'\t%4.4f',PTanova(i,j,t,1));
            end
        end
    end
end
    
fclose(fp);





function SaveFigure(obj, event,Device,mask,comment)

if (isempty(mask))
    print(Device);
else

    [filename, pathname] = uiputfile(mask,comment);
    if isequal(filename,0) || isequal(pathname,0)
        return
    end
    print(Device,fullfile(pathname, filename));
end



function ShowTanovaBtCallbackKid(obj, event, in1,in2)
ph = get(obj,'Parent');
if (ph ~= 0)
    ShowTanovaBtCallback(ph, event, in1,in2)
end


function ShowTanovaBtCallback(obj, event, in1,in2)

pos = get(obj,'CurrentPoint');


%fh = figure(103);
fh = gcf;

tanh = findobj(obj,'Tag','TanovaP');
x = get(tanh,'XData');
if isempty(x)
    x = 1;
end
delta = abs(x - pos(1,1));

[mn,idx] = min(delta);
d = get(get(obj,'Parent'),'UserData');

if in1 == 0
    idx = 1;
end

ch = findobj(fh,'Tag','Cursor');
for i = 1:numel(ch)
    set(ch(i),'XData',[x(idx) x(idx)]);
end

sh = findobj(fh,'Tag','subplot');

for i = 1:numel(sh)
    set(get(sh(i),'Title'),'FontWeight','light');
end

set(get(obj,'Title'),'FontWeight','bold');

if (d.Averaged == 0)
    txt = cell(numel(sh),1);
    for i = 1:numel(sh)
        tanh = findobj(sh(i),'Tag','TanovaP');
        y = get(tanh,'YData');
        txt{i} = sprintf('%s: %3.3f (%1.1f %s)',get(sh(i),'UserData'),y(idx),x(idx),d.axislabel);
        set(get(sh(i),'Title'),'String',txt{i},'Interpreter','none','FontSize',8);
    end
else
    txt = cell(numel(sh),1);
    for i = 1:numel(sh)
        tanh = findobj(sh(i),'Tag','TanovaP');
        y = get(tanh,'YData');
        txt{i} = sprintf('%s: %3.3f',get(sh(i),'UserData'),y);
        set(get(sh(i),'Title'),'String',txt{i},'Interpreter','none','FontSize',8);
    end
    
end


if in1 == 1
    Mean2Show = d.gmeans(:,:,:,idx);
    GFP2Show  = d.gmeansGFP(:,:,:,idx);
   
else
    Mean2Show = d.means(:,:,:,idx);
    GFP2Show  = d.meansGFP(:,:,:,idx);
end

if in2 == 1
    Mean2Show = mean(Mean2Show,2);
    GFP2Show = mean(GFP2Show,2);
end

if in2 == 2
    if d.ContF1 == false
        idx_map = unique(d.SelDesign(:,1));
    
        nm    = zeros(size(Mean2Show,1),numel(idx_map),size(Mean2Show,3));
        nmGFP = zeros(size(GFP2Show,1),numel(idx_map),size(GFP2Show,3));
        for i = 1:numel(idx_map)
            nm(:,i,:)    = squeeze(mean(Mean2Show(:,d.SelDesign(:,1) == idx_map(i),:),2));
            nmGFP(:,i,:) = squeeze(mean(GFP2Show(:,d.SelDesign(:,1)  == idx_map(i),:),2));
        end
        
        Mean2Show = nm;
        GFP2Show = nmGFP;
    else
        if in1 == 2
            Mean2Show = d.cmeans(:,:,:,idx);
            GFP2Show  = d.cmeansGFP(:,:,:,idx);
        end
    end
end

if in2 == 3
    idx_map = unique(d.SelDesign(:,2));
    nm = zeros(size(Mean2Show,1),numel(idx_map),size(Mean2Show,3));
    nmGFP = zeros(size(GFP2Show,1),numel(idx_map),size(GFP2Show,3));
    for i = 1:numel(idx_map)
        nm(:,i,:)    = squeeze(mean(Mean2Show(:,d.SelDesign(:,2) == idx_map(i),:),2));
        nmGFP(:,i,:) = squeeze(mean(GFP2Show(:,d.SelDesign(:,2) == idx_map(i),:),2));
    end
    Mean2Show = nm;
    GFP2Show = nmGFP;
end

if in2 == 4
    nLevel1 = max(d.SelDesign(:,1));
    mn = min(d.SelDesign(:,2));
    iDesign = d.SelDesign(:,1) + (nLevel1+1) * (d.SelDesign(:,2)-mn +1);
    idx_map = unique(iDesign);
    % fix this
    nm = zeros(size(Mean2Show,1),numel(idx_map),size(Mean2Show,3));
    nmGFP = zeros(size(GFP2Show,1),numel(idx_map),size(GFP2Show,3));
    for i = 1:numel(idx_map)
        nm(:,i,:)    = squeeze(mean(Mean2Show(:,iDesign == idx_map(i),:),2));
        nmGFP(:,i,:) = squeeze(mean(GFP2Show(:,iDesign == idx_map(i),:),2));
    end
    Mean2Show = nm;
    GFP2Show = nmGFP;
end

ng = size(Mean2Show,1);
nc = size(Mean2Show,2);

d.lbl = cell(nc,1);

switch in2
    case 1
        d.lbl = {'Main effect'};
    case 2
        if d.ContF1
            d.lbl{1} = [d.strF1 '+'];
            d.lbl{2} = [d.strF1 '-'];
        else
           [d.lbl{1:nc}] = deal(d.DLabels1.Label);
        end
    case 3
        [d.lbl{1:nc}] = deal(d.DLabels2.Label);
    case 4
        for i = 1:numel(idx_map)
            l1 = rem(idx_map(i),nLevel1+1);
            l2 = (idx_map(i) - l1) / (nLevel1+1) +mn -1;
            [ll1c{1:numel(d.DLabels1)}] = deal(d.DLabels1.Level);
            ll1 = cell2mat(ll1c);
            idx1 = find(ll1 == l1,1);
            [ll2c{1:numel(d.DLabels2)}] = deal(d.DLabels2.Level);
            ll2 = cell2mat(ll2c);
            idx2 = find(ll2 == l2,1);
            d.lbl{i} = sprintf('%s\n%s',d.DLabels1(idx1).Label,d.DLabels2(idx2).Label);
            
        end
%         lbl1 = cell(numel(d.DLabels1),1);
%         lbl2 = cell(numel(d.DLabels1),1);
%         [i1,j1,k1] = unique(d.SelDesign(:,1));
%         [dummy,j2,k2] = unique(d.SelDesign(:,2));
%         [lbl1{1:numel(d.DLabels1)}] = deal(d.DLabels1.Label);
%         [lbl2{1:numel(d.DLabels2)}] = deal(d.DLabels2.Label);
%         for c = 1:nc
%             d.lbl{c} = sprintf('%s\n%s',lbl1{k1(c)},lbl2{k2(c)});
%         end
end


mapScale = max(abs(Mean2Show(:))) / 8;
cntlbl2 = {'+','-'};

if isfield(d,'Channel') && d.DoGFP == 0
%    disp('Showing mean maps')
    for g = 1:ng
        for c = 1:nc
            mh = subplot(2*ng,2*nc,(g-1)*2*nc+nc +c);
            set(mh,'ActivePositionProperty','outerposition','Box','on');
%            mh = subaxis(2*ng,2*nc,(g-1)*2*nc+nc +c,'SpacingHoriz', 0.03,'SpacingVert', 0.03, 'Padding', 0, 'Margin', 0);
            axis equal
            dt = squeeze(Mean2Show(g,c,:));
            dt = dt(:);

            pos = get(mh,'Position');
            opos = get(mh,'OuterPosition');
            set(mh,'OuterPosition',[pos(1) pos(2) pos(3) opos(4)]);

            if ng > 1
                if d.ContBetween == false
                    tit = (sprintf('%s\n%s',char(d.GLabels{g}),char(d.lbl{c})));
                    stit = (sprintf('%s - %s',char(d.GLabels{g}),char(d.lbl{c})));
                else
                    tit = (sprintf('%s%s\n%s',char(d.IndName),char(cntlbl2{g}),char(d.lbl{c})));
                    stit = (sprintf('%s%s - %s',char(d.IndName),char(cntlbl2{g}),char(d.lbl{c})));
                end
            else
                if d.ContF1 == false
                    tit = (sprintf('%s\n%s',d.MainGroupLabel,char(d.lbl{c})));
                    stit = (sprintf('%s - %s',d.MainGroupLabel,char(d.lbl{c})));
                else
                    tit = (sprintf('%s',char(d.lbl{c})));
                    stit = (sprintf('%s',char(d.lbl{c})));
                end
            end
           
            
            RaguDSPMap(dt,d.Channel,d.MapStyle,mapScale,'Title',stit); % ,'NoScale','Resolution',3);
            th = title(tit);
            set(th,'Interpreter','none','FontSize',8);
        end
    end
end

if d.DoGFP == 0

    maps = reshape(Mean2Show,size(Mean2Show,1)*size(Mean2Show,2),size(Mean2Show,3));
    maps = maps - repmat(mean(maps),size(maps,1),1);
    cov = maps' * maps;

    opt.disp = 0;
    [v,de] = eigs(cov,2,'LM',opt);
    d.pro = maps*v;
    
%    ops = statset('Display','off','TolFun',0.01);
% %   squareform(pdist(maps,'euclidean'))
    TwoDimsToMap = size(maps,1) > 2;
%    if TwoDimsToMap == true
%        d.pro = mdscale(pdist(maps,'euclidean'),2,'Options',ops,'Criterion','strain');
%    else
%        d.pro = mdscale(pdist(maps,'euclidean'),1,'Options',ops,'Criterion','strain');
%        d.pro(:,2) = 0;
%    end
    
    v = NormDimL2(maps'*d.pro,1);

    d.pro = reshape(d.pro,size(Mean2Show,1),size(Mean2Show,2),2);
    d.timePT = x(idx);
    d.in1 = in1;

    set(get(obj,'Parent'),'UserData',d);

    sh = ShowMDS(d,false);

    spos = get(sh,'Position');

    col = lines();
    psymb = '.ox+*sdv^<>ph';
    subplot('Position',[0.88 0.15 0.1 0.35]);
    axis off
    cla
    if d.ContBetween == false || in1 == 1
        if ng > 1 && in1 > 1
%            patch([0 8 8 0],[20 20 (19 - numel(d.GLabels)) (19 - numel(d.GLabels))],[1 1 1]);
            hold on
            for g = 1:numel(d.GLabels)
                plot(1,20-g,psymb(g),'Color',col(g,:));
                text(2,20-g,d.GLabels{g},'Interpreter','none','VerticalAlignment','middle','FontSize',8,'Color',col(g,:));
            end
            hold off
            axis([0 8 0 20]);
            axis off
            hold off
        end
    else
 %       subplot('Position',[0.88 0.15 0.1 0.35]);
        linecolors = lines(size(d.pro,2));
%        cla
%        patch([0 8 8 0],[20 20 (19 - numel(d.lbl)) (19 - numel(d.lbl))],[1 1 1]);
        hold on
  
        for i = 1:numel(d.lbl)
            plot(1,20-i,'+','Color',linecolors(i,:));
            text(2,20-i,d.lbl{i},'Interpreter','none','VerticalAlignment','middle','FontSize',8,'Color',linecolors(i,:));
           hold on
        end
        axis([0 8 0 20]);
        axis off
    end

    if isfield(d,'Channel')
        dlt = 0.03;
%        disp('Showing MDS maps');
        subplot('Position',[(spos(1)        ) (spos(2)-3*dlt) 2*dlt 2*dlt]);
        RaguDSPMap(-v(:,1),d.Channel,d.MapStyle); %,'NoScale','Resolution',3);
    
        subplot('Position',[(spos(1)+spos(3)-2*dlt) (spos(2)-3*dlt) 2*dlt 2*dlt]);
        RaguDSPMap( v(:,1),d.Channel,d.MapStyle);% ,'NoScale','Resolution',3);

        if TwoDimsToMap == true
            subplot('Position',[(spos(1)-2*dlt) (spos(2)) 2*dlt 2*dlt]);
            RaguDSPMap(-v(:,2),d.Channel,d.MapStyle); % ,'NoScale','Resolution',3);
    
            subplot('Position',[(spos(1)-2*dlt) (spos(2)+spos(4)-2*dlt) 2*dlt 2*dlt]);
            RaguDSPMap( v(:,2),d.Channel,d.MapStyle); % ,'NoScale','Resolution',3);
        else
            subplot('Position',[(spos(1)-2*dlt) (spos(2)) 2*dlt 2*dlt]);
            cla
            axis off
            subplot('Position',[(spos(1)-2*dlt) (spos(2)+spos(4)-2*dlt) 2*dlt 2*dlt]);
            cla
            axis off
        end
    end
else
    subplot(224);
%    gfp = squeeze(std(Mean2Show,1,3));
    gfp = GFP2Show;
    bh = bar(gfp);
    bcm = lines(32);
    xlabel('Condition');
    ylabel('GFP');

    for i = 1:numel(bh)
        set(bh(i),'FaceColor',bcm(i,:));
    end
    h = gca;
    
    if size(gfp,1) > 1
        if d.ContBetween == 0
            set(h,'XTickLabel',d.GLabels,'TickLabelInterpreter','none');
        else
            set(h,'XTickLabel',{[d.IndName,'+'],[d.IndName,'-']},'TickLabelInterpreter','none');
        end
        if (size(gfp,2) > 1)
            legend(d.lbl,'Location','NorthOutside','interpreter','none');
            
        end
    else
        set(h,'XTickLabel',d.lbl);
        if in2 == 4
            rotateticklabel(h,90);
        end
    end
end


function handle = ShowMDS(d,newFig)

if nargin < 2
    newFig = false;
end
if newFig == true
    handle = axes;
else
    handle = subplot('Position',[0.60 0.15 0.25 0.35]);
end

cntlbl = {'+','-'};
col = lines();
hold off

psymb = '.ox+*sdv^<>ph';

if d.ContBetween == false || d.in1 == 1
    for g = 1:size(d.pro,1)
        plot(squeeze(d.pro(g,:,1)),squeeze(d.pro(g,:,2)),psymb(g),'Color',col(g,:));
        if (d.ContBetween == false)
            axis([-d.mx *d.mdsScaleFactor d.mx * d.mdsScaleFactor -d.mx * d.mdsScaleFactor d.mx * d.mdsScaleFactor]);
        else
            axis([-d.mxgr * d.mdsScaleFactor d.mxgr * d.mdsScaleFactor -d.mxgr * d.mdsScaleFactor d.mxgr * d.mdsScaleFactor]);
        end
        hold on
        th = text(squeeze(d.pro(g,:,1)),squeeze(d.pro(g,:,2)),d.lbl,'Color',col(g,:));
        set(th,'Interpreter','none','VerticalAlignment','bottom','FontSize',8);
    end

else
    pro_x = d.pro(:,:,1);
    pro_y = d.pro(:,:,2);

    plot(pro_x,pro_y,'-');
    axis equal
    for i = 1:size(pro_y,2)
       text(pro_x(:,i),pro_y(:,i),cntlbl,'Interpreter','none','VerticalAlignment','bottom','FontSize',8,'Color',col(i,:));
    end
        
    axis([-d.mx*d.mdsScaleFactor d.mx*d.mdsScaleFactor -d.mx*d.mdsScaleFactor d.mx*d.mdsScaleFactor]);
end

if (d.Averaged == 0)
    title(sprintf('MDS at %1.1f %s',d.timePT,d.axislabel));
else
    title('MDS');
end


function RescaleMDS(obj,eventdata,fact,newFig)
d = get(get(obj,'Parent'),'UserData');
d.mdsScaleFactor = d.mdsScaleFactor * fact;
set(gcf,'UserData',d);
if newFig == true
    figure
end
ShowMDS(d,newFig);


% function ExportAMap(obj,eventdata)
% 
% ud = get(obj,'UserData');
% name = ud(1).Name';
% name = strtrim(name(:)');
% [fn, pn] = uiputfile([name '.txt'], 'Save map to textfile');
% 
% if fn == 0 
%     return
% end
% 
% fp = fopen([pn fn],'wt');
% for i = 1:numel(ud.Channel)
%     fprintf(fp,'\t%s',ud.Channel(i).Name);
% end
% fprintf(fp,'\n%s',ud.Name);
% 
% for i = 1:numel(ud.Channel)
%     fprintf(fp,'\t%f',ud.Map(i));
% end

% function PlotAMap(obj,eventdata)
% 
% ud = get(obj,'UserData');
% subplot('Position',[0.1 0.2 0.8 0.8]);
% %for i = 1:numel(montage)
% %    lbl{i} = montage(i).Name;
% %end
% title(ud.Name);
% RaguDSPMap(ud.Map,ud.Channel,ud.Style,'NoScale','Resolution',1,'Step',ud.Scale);%,'Label',lbl');
% 
% subplot('Position',[0.3 0.1 0.4 0.05]);
% 
% dspCMapColorbar(ud.Scale,'br');
% 


