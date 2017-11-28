function varargout = Ragu_Cormat(varargin)
% Ragu_Cormat M-file for Ragu_Cormat.fig
%      Ragu_Cormat, by itself, creates a new Ragu_Cormat or raises the existing
%      singleton*.
%
%      H = Ragu_Cormat returns the handle to a new Ragu_Cormat or the handle to
%      the existing singleton*.
%
%      Ragu_Cormat('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in Ragu_Cormat.M with the given input arguments.
%
%      Ragu_Cormat('Property','Value',...) creates a new Ragu_Cormat or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Ragu_Cormat_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Ragu_Cormat_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Ragu_Cormat

% Last Modified by GUIDE v2.5 07-Jan-2015 20:31:51

% Begin initialization code - DO NOT EDIT

global CURRENTVIEW
CURRENTVIEW = 'Cormat';

gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Ragu_Cormat_OpeningFcn, ...
                   'gui_OutputFcn',  @Ragu_Cormat_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before Ragu_Cormat is made visible.
function Ragu_Cormat_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Ragu_Cormat (see VARARGIN)

% Choose default command line output for Ragu_Cormat
handles.output = hObject;
handles.CurrentView = 'Cormat';

ud.CondLabels  = varargin{4}.conds;
ud.ContBetween = varargin{4}.ContBetween;

if ud.ContBetween == 0
    ud.GroupLabels = varargin{4}.GroupLabels;
else
    ud.GroupLabels = {'All subjects'};
end

set(handles.GroupPosList,'String',ud.GroupLabels,'Min',0,'Max',100);
set(handles.GroupNegList,'String',ud.GroupLabels,'Min',0,'Max',100);
set(handles.CondPosList,'String',ud.CondLabels,'Min',0,'Max',100);
set(handles.CondNegList,'String',ud.CondLabels,'Min',0,'Max',100);

ud.GPosIdx = 1;
ud.GNegIdx = 1;
ud.CPosIdx = 1;
ud.CNegIdx = 1;

ud.Montage = varargin{4}.Channel;
ud.TimeOnset = varargin{4}.TimeOnset;
ud.DeltaX = varargin{4}.DeltaX;
ud.V = varargin{4}.V;
ud.IndFeature = varargin{4}.IndFeature;
ud.MapStyle = varargin{4}.MapStyle;
ud.txtX = varargin{4}.txtX;

ud.time = (0:(size(ud.V,4)-1));
ud.time = ud.time * ud.DeltaX + ud.TimeOnset;

axis(handles.axes_Cormat,[1,size(ud.V,4),1,size(ud.V,4)]); 

ytick = get(handles.axes_Cormat,'YTick');
xtick = get(handles.axes_Cormat,'XTick');
set(handles.axes_Cormat,'YTickLabel',ytick * ud.DeltaX + ud.TimeOnset,'XTickLabel',xtick * ud.DeltaX + ud.TimeOnset);

set(handles.Ragu_Cormat,'UserData',ud);

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes Ragu_Cormat wait for user response (see UIRESUME)
uiwait(handles.Ragu_Cormat);

function TMapCellSelection(hObject, Indices, handles,tableID)

ud = get(handles.Ragu_Cormat,'UserData');

%set(handles.sLoreta,'Enable','off');

switch tableID
    case 1
%        ud.CPosIdx = eventdata.Indices(:,1);
        ud.CPosIdx = Indices(:,1);
    case 2  
%        ud.CNegIdx = eventdata.Indices(:,1);
        ud.CNegIdx = Indices(:,1);
    case 3  
%        ud.GPosIdx = eventdata.Indices(:,1);
        ud.GPosIdx = Indices(:,1);
    case 4  
%        ud.GNegIdx = eventdata.Indices(:,1);
        ud.GNegIdx = Indices(:,1);
end

set(handles.Ragu_Cormat,'UserData',ud);


function res = GetSubjectSelection(ud,idx)
res = zeros(size(ud.V,1),1);

if ud.ContBetween == 0
    for i = 1:numel(idx)
        res(ud.IndFeature == idx(i)) = 1;
    end
else
    res = ~isnan(ud.IndFeature);
end

function res = PreprocessedMean(d)

res = squeeze(mean(mean(d,4),2));

function res = IsIndexDifferent(i1,i2)

if (numel(i1) ~= numel(i2))
    res = true;
    return;
end

i1s = sort(i1(:));
i2s = sort(i2(:));

if sum(i1s ~= i2s) == 0
    res = false;
    return;
end
res = true;

function DoTheCormat(handles)

ud = get(handles.Ragu_Cormat,'UserData');

if isempty(ud.CPosIdx) + isempty(ud.GPosIdx) + isempty(ud.CNegIdx) + isempty(ud.GNegIdx) > 0
    return
end

S1ToUse = GetSubjectSelection(ud,ud.GPosIdx);
S2ToUse = GetSubjectSelection(ud,ud.GNegIdx);
    
Pos = squeeze(mean(mean(ud.V(S1ToUse == 1,ud.CPosIdx,:,:),1),2));
Neg = squeeze(mean(mean(ud.V(S2ToUse == 1,ud.CNegIdx,:,:),1),2));

cm = corr(Neg,Pos);
axes(handles.axes_Cormat);

h = imagesc(cm',[-1 1]);
axis xy
colormap(BlueRed);
ytick = get(handles.axes_Cormat,'YTick');
xtick = get(handles.axes_Cormat,'XTick');

axud.Pos = Pos;
axud.Neg = Neg;
axud.Mnt = ud.Montage;
axud.Style = ud.MapStyle;
axud.DeltaX = ud.DeltaX;
axud.off = ud.TimeOnset;
axud.txtX = ud.txtX;
set(handles.axes_Cormat,'YTickLabel',ytick * ud.DeltaX + ud.TimeOnset,'XTickLabel',xtick * ud.DeltaX + ud.TimeOnset);
xlabel(handles.axes_Cormat,ud.txtX);
ylabel(handles.axes_Cormat,ud.txtX);
title(handles.axes_Cormat,'Correlation matrix');
set(h,'ButtonDownFcn',{@Cormat_ButtonDownFcn,handles,axud});

% --- Outputs from this function are returned to the command line.
function varargout = Ragu_Cormat_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
if isempty(handles)
    varargout{1} = [];
else
    ud = get(handles.Ragu_Cormat,'UserData');
    set(handles.output,'UserData',ud);
    varargout{1} = handles.output;
end


function MapBtCallback(obj, event, data,montage,style,step)

figure;
subplot('Position',[0.1 0.2 0.8 0.8]);
%for i = 1:numel(montage)
%    lbl{i} = montage(i).Name;
%end
RaguDSPMap(data,montage,style,'NoScale','Resolution',1,'Step',step);%,'Label',lbl');
subplot('Position',[0.3 0.1 0.4 0.05]);
dspCMapColorbar(step,'br');


% --- Executes on button press in ButtonDone.
function ButtonDone_Callback(hObject, eventdata, handles)

uiresume(handles.Ragu_Cormat);

% --- Executes on selection change in GroupPosList.
function GroupPosListCM_Callback(hObject, eventdata, handles)


%eventdata.Indices = get(hObject,'Value')';
%TMapCellSelection(hObject, eventdata, handles,3);
Indices = get(hObject,'Value')';
TMapCellSelection(hObject, Indices, handles,3);

% --- Executes during object creation, after setting all properties.
function GroupPosListCM_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in GroupNegList.
function GroupNegListCM_Callback(hObject, eventdata, handles)

%eventdata.Indices = get(hObject,'Value')';
%TMapCellSelection(hObject, eventdata, handles,4);
Indices = get(hObject,'Value')';
TMapCellSelection(hObject, Indices, handles,4);



% --- Executes during object creation, after setting all properties.
function GroupNegListCM_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in CondPosList.
function CondPosListCM_Callback(hObject, eventdata, handles)

%eventdata.Indices = get(hObject,'Value')';
%TMapCellSelection(hObject, eventdata, handles,1);
Indices = get(hObject,'Value')';
TMapCellSelection(hObject, Indices, handles,1);

% --- Executes during object creation, after setting all properties.
function CondPosListCM_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in CondNegList.
function CondNegListCM_Callback(hObject, eventdata, handles)

%eventdata.Indices = get(hObject,'Value')';
%TMapCellSelection(hObject, eventdata, handles,2);
Indices = get(hObject,'Value')';
TMapCellSelection(hObject, Indices, handles,2);


% --- Executes during object creation, after setting all properties.
function CondNegListCM_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function editTitle_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function editTitle_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in ButtonShow.
function ButtonShow_Callback(hObject, eventdata, handles)
% hObject    handle to ButtonShow (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
DoTheCormat(handles);


% --- Executes on mouse press over axes background.
function axes_Cormat_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to axes_Cormat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


function Cormat_ButtonDownFcn(hObject, eventdata,handles,ud)
par = get(hObject,'Parent');
pt = round(get(par,'CurrentPoint'));
nf = size(ud.Pos,2);

pt(pt > nf) = nf;
pt(pt <  1) = 1;

c = corr(ud.Pos(:,pt(1,2)),ud.Neg(:,pt(1,1)));

axes(handles.axesMap1);
RaguDSPMap(ud.Pos(:,pt(1,2)),ud.Mnt,ud.Style,'NoScale','Resolution',1);
axes(handles.axesMap2);
RaguDSPMap(ud.Neg(:,pt(1,1)),ud.Mnt,ud.Style,'NoScale','Resolution',1);


set(handles.textTimeY,'String',sprintf('%2.2f%s',pt(1,2) * ud.DeltaX + ud.off,ud.txtX));
set(handles.textTimeX,'String',sprintf('%2.2f%s',pt(1,1) * ud.DeltaX + ud.off,ud.txtX));
title(handles.axes_Cormat,['Correlation: ' sprintf('%3.3f',c)]);

% --- Executes during object creation, after setting all properties.
function Ragu_Cormat_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Ragu_Cormat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes when user attempts to close Ragu_Cormat.
function Ragu_Cormat_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to Ragu_Cormat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
delete(hObject);


% --- Executes during object deletion, before destroying properties.
function Ragu_Cormat_DeleteFcn(hObject, eventdata, handles)
% hObject    handle to Ragu_Cormat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes when Ragu_Cormat is resized.
function Ragu_Cormat_ResizeFcn(hObject, eventdata, handles)
% hObject    handle to Ragu_Cormat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over ButtonShow.
function ButtonShow_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to ButtonShow (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
