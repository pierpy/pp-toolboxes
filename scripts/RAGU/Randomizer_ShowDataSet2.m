function Randomizer_ShowDataSet2()

% Copyright 2009-2011 Thomas Koenig
% distributed under the terms of the GNU AFFERO General Public License

uitableNames = uitable;
out = get(gcf,'UserData');
pos = [0.05 0.10 0.45 0.85];
set(uitableNames,'Data',out.Names,'ColumnName',out.conds,'Units','Normalized','Position',pos,'BackgroundColor',[1 1 1]);
set(uitableNames,'Tag','DataTable');
set(uitableNames,'Units','points');

%p = floor((get(uitableNames,'Position') -55)/ numel(out.conds));
%pc = repmat({p(3)},1,numel(out.conds));
%pc = repmat({'auto'},1,numel(out.conds))
%set(uitableNames,'ColumnWidth','auto');
set(uitableNames,'Units','Normalized');
ext = get(uitableNames,'Extent');

pos(3) = min([ext(3),pos(3)]);
pos(4) = min([ext(4),pos(4)]);
pos(2) = 0.95 - pos(4);

handles = guidata(gcf);
handles.CurrentView = 'Data';
guidata(gcf,handles);

uicontrol('Style', 'text', 'String', 'Select datasets to show','Units','Normalized','Position', [0.05 0.96 0.45 0.03],'BackgroundColor',[1 1 1],'Tag','DataTable');

%set(uitableNames,'Position',pos);


d.time = (0:(size(out.V,4)-1)) * out.DeltaX + out.TimeOnset;
d.txtX = out.txtX;
d.MapStyle = out.MapStyle;
d.pos = d.time(1);
d.max = max(abs(out.V(:)));
d.fact = 1;
hPlotArea = subplot('Position',[0.55 0.4 0.4 0.55]);
d.hMapArea  = subplot('Position',[0.68 0.06 0.2 0.2]);
axis(hPlotArea,'off');
axis(d.hMapArea,'off');
set(uitableNames,'UserData',d);
%UpdateButtonHandle = uicontrol('Style', 'pushbutton', 'String', 'Show','Units','Normalized','Position', [0.05 pos(2)-0.08 pos(3) 0.075], 'Callback', @ShowTheData,'Tag','DataTable');
UpdateButtonHandle = uicontrol('Style', 'pushbutton', 'String', 'Show','Units','Normalized','Position', [0.05 0.05 0.45 0.05], 'Callback', @ShowTheData,'Tag','DataTable');
PlusButtonHandle = uicontrol('Style', 'pushbutton', 'String', '+','Units','Normalized','Position', [0.955 0.90 0.04 0.04], 'Callback', @ButtonPlus_Callback,'Tag','DataTable');
MinusButtonHandle = uicontrol('Style', 'pushbutton', 'String', '-','Units','Normalized','Position', [0.955 0.85 0.04 0.04], 'Callback', @ButtonMinus_Callback,'Tag','DataTable');


eventdata.Indices = [1 1];
handles.uitableNames = uitableNames;
handles.PlotArea = hPlotArea;
handles.MapArea  = d.hMapArea;
handles.Indices= eventdata.Indices;
handles.ShowButton = UpdateButtonHandle;
set(handles.PlotArea,'UserData',d);
ButtonUD.Indices = eventdata.Indices;
ButtonUD.PlotArea = hPlotArea;
set(uitableNames,'UserData',handles);
set(UpdateButtonHandle,'UserData',ButtonUD);
set(PlusButtonHandle,'UserData',ButtonUD);
set(MinusButtonHandle,'UserData',ButtonUD);

uitableNames_CellSelectionCallback(uitableNames,eventdata,handles);
set(uitableNames,'CellSelectionCallback',@uitableNames_CellSelectionCallback);


% --- Executes when selected cell(s) is changed in uitableNames.
function uitableNames_CellSelectionCallback(hObject,eventdata,handles)
% hObject    handle to uitableNames (see GCBO)
% eventdata  structure with the following fields (see UITABLE)
%	Indices: row and column indices of the cell(s) currently selecteds
% handles    structure with handles and user data (see GUIDATA)
handles = get(hObject,'UserData');
ud = get(handles.ShowButton,'UserData');
ud.Indices = eventdata.Indices;
set(handles.ShowButton,'UserData',ud);
cla(handles.PlotArea);
axis(handles.MapArea);
title('');
cla(handles.MapArea);
axis(handles.PlotArea,'off');
% handles.Indices = eventdata.Indices;

function ShowTheData(hObject, eventdata)
out = get(gcf,'UserData');
ud = get(hObject,'UserData');
d = get(ud.PlotArea,'UserData');
axis(ud.PlotArea,'on');

%eventdata
%eventdata.Indices = ud.Indices;
%d = get(ud.PlotArea,'UserData');

d.MMaps = zeros(size(out.V,3),size(out.V,4));
if (isfield(out,'Channel'))
    d.Channel = out.Channel;
end

for i = 1:size(ud.Indices,1)
    d.MMaps = d.MMaps + squeeze(out.V(ud.Indices(i,1),ud.Indices(i,2),:,:));
end
d.MMaps = d.MMaps / size(ud.Indices,1);
set(ud.PlotArea,'UserData',d);
%ud.Indices = eventdata.Indices;
set(hObject,'UserData',ud);
ShowCurves(ud.PlotArea);


function ShowCurves(hObject)

d = get(hObject,'UserData');
axes(hObject);

plot(d.time,d.MMaps,'-k');
hold on
d.curh = plot([d.pos(1) d.pos(1)],[-10000 10000],'-r');
hold off
if (numel(d.time) > 1)
    axis([min(d.time) max(d.time) (-d.max*1.05) (d.max*1.05)]);
end
title('Mean of selected cells');

set(hObject,'UserData',d);

kids = get(hObject,'Children');

for c = 1:numel(kids)
    set(kids(c),'ButtonDownFcn',@ShowMapKid);
end

set(hObject,'UserData',d);
set(hObject,'ButtonDownFcn',@ShowMap);
%subplot('Position',[0.68 0.055 0.2 0.201]);
%axis off

%ShowMap(hObject,[d.pos d.pos]);



function ShowMapKid(hObject,eventdata,handles)
ShowMap(get(hObject,'Parent'));



function ShowMap(hObject,pos)


if (nargin < 2) || isempty(pos)
    pos = get(hObject,'CurrentPoint');
else
    if verLessThan('8.0','matlab') == false
        pos = pos.IntersectionPoint(1:2);
    end
end

d = get(hObject,'UserData');

delta = abs(d.time - pos(1,1));
[mn,idx] = min(delta);

d.pos = d.time(idx);
set(d.curh,'XData',[d.pos d.pos]);
set(hObject,'UserData',d);
if (isfield(d,'Channel'))
    axes(d.hMapArea);
%    subplot('Position',[0.68 0.06 0.2 0.2]);
    RaguDSPMap(d.MMaps(:,idx),d.Channel,d.MapStyle,'NoScale','Resolution',3);
    title(sprintf('%3.1f %s',d.pos,d.txtX),'FontSize',10);
end



% --- Executes on button press in ButtonPlus.
function ButtonPlus_Callback(hObject, eventdata, handles)
% hObject    handle to ButtonPlus (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ud = get(hObject,'UserData');
d = get(ud.PlotArea,'UserData');
d.max = d.max / 1.5;
set(ud.PlotArea,'UserData',d);

ShowCurves(ud.PlotArea);


% --- Executes on button press in ButtonMinus.
function ButtonMinus_Callback(hObject,eventdata)
% hObject    handle to ButtonMinus (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ud = get(hObject,'UserData');
d = get(ud.PlotArea,'UserData');
d.max = d.max * 1.5;
set(ud.PlotArea,'UserData',d);

ShowCurves(ud.PlotArea);


