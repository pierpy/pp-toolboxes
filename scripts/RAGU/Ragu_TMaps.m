function varargout = Ragu_TMaps(varargin)
% RAGU_TMAPS M-file for Ragu_TMaps.fig
%      RAGU_TMAPS, by itself, creates a new RAGU_TMAPS or raises the existing
%      singleton*.
%
%      H = RAGU_TMAPS returns the handle to a new RAGU_TMAPS or the handle to
%      the existing singleton*.
%
%      RAGU_TMAPS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in RAGU_TMAPS.M with the given input arguments.
%
%      RAGU_TMAPS('Property','Value',...) creates a new RAGU_TMAPS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Ragu_TMaps_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Ragu_TMaps_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Ragu_TMaps

% Last Modified by GUIDE v2.5 23-Jun-2011 18:00:38

% Begin initialization code - DO NOT EDIT
gui_Singleton = 0;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Ragu_TMaps_OpeningFcn, ...
                   'gui_OutputFcn',  @Ragu_TMaps_OutputFcn, ...
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


% --- Executes just before Ragu_TMaps is made visible.
function Ragu_TMaps_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Ragu_TMaps (see VARARGIN)

% Choose default command line output for Ragu_TMaps
handles.output = hObject;

ud = varargin{4};

subplot(handles.MapPos);
RaguDSPMap(ud.MPos,ud.Montage,2,'NoScale','Resolution',2,'Step',ud.MRange,'tValue',ud.tPos,'Title','Mean positive');
%set(get(handles.MapPos,'Children'),'ButtonDownFcn',{@MapBtCallback,ud.MPos,ud.Montage,ud.MapStyle,ud.MRange});
subplot(handles.Meancolorbar1);
dspCMapColorbar(ud.MRange,'br');

subplot(handles.MapNeg);
RaguDSPMap(ud.MNeg,ud.Montage,2,'NoScale','Resolution',2,'Step',ud.MRange,'tValue',ud.tNeg,'Title','Mean negative');
%set(get(handles.MapNeg,'Children'),'ButtonDownFcn',{@MapBtCallback,ud.MNeg,ud.Montage,ud.MapStyle,ud.MRange});
subplot(handles.Meancolorbar2);
dspCMapColorbar(ud.MRange,'br');

if ~isempty(ud.t)
    subplot(handles.TMap);
    if ud.ShowPValues == false
        set(handles.text3,'String','t-map');
        RaguDSPMap(ud.t,ud.Montage,2,'NoScale','Resolution',1,'Step',ud.TRange,'Title','t-values');
%        set(get(handles.TMap,'Children'),'ButtonDownFcn',{@MapBtCallback,ud.t,ud.Montage,ud.MapStyle,ud.TRange});
        subplot(handles.Tcolorbar);
        cp = get(handles.Tcolorbar,'Position');
        set(handles.Tcolorbar,'Position',[cp(1),22,cp(3),3]);
        dspCMapColorbar(ud.TRange,'br');
    
    else
        LevelListNeg = [0.001 0.005 0.01 0.02 0.05 0.1];
        LevelListPos = [0.1 0.05 0.02 0.01 0.005 0.001];
        
%        tCumDist = [(LevelListNeg / 2) 0 (1 -LevelListPos / 2)];
%        LevelList = tinv(tCumDist,ud.df);

        LevelList= [tinv(LevelListNeg / 2,ud.df) 0 tinv(1 -LevelListPos / 2,ud.df)];

        [maxT,maxT_idx] = max(ud.t);
        [minT,minT_idx] = min(ud.t);
        
        lbl = {};
        lblidx = [];
        if maxT > tinv(0.975,ud.df)
            lblidx(1) = maxT_idx;
            lbl{numel(lblidx)} = '+';
        end
        
        if minT < tinv(0.025,ud.df)
            lblidx = [lblidx, minT_idx];
            lbl{numel(lblidx)} = '-';
        end
        RaguDSPMap(ud.t,ud.Montage,2,'NoScale','Resolution',1,'LevelList',LevelList,'Label',lbl,'LabelIndex',lblidx,'Title','p-value');
%        set(get(handles.TMap,'Children'),'ButtonDownFcn',{@PMapBtCallback,ud.t,ud.Montage,LevelList,lbl,lblidx,[-LevelListNeg 0 LevelListPos]});
        set(handles.text3,'String','p-map');
        subplot(handles.Tcolorbar);
        cp = get(handles.Tcolorbar,'Position');
        set(handles.Tcolorbar,'Position',[cp(1),22,cp(3),3]);
        dspPMapColorbar([-LevelListNeg 0 LevelListPos]);
    end
else
    set(handles.text3,'String','no contrast');
end

set(handles.editPos,'String',ud.PosString);
set(handles.editNeg,'String',ud.NegString);

set(handles.figure1,'Name',ud.Title);

if isempty(ud.t)
    res =  [ud.MPos;ud.MNeg];
else
    res =  [ud.MPos;ud.MNeg;ud.t];
end

CLabel = cell(numel(ud.Montage),1);


for i = 1:numel(ud.Montage)
    if isfield(ud.Montage(i),'Name')
        CLabel{i} = ud.Montage(i).Name;
    else
        CLabel{i} = sprintf('Ch%i',i);
    end
end
    
set(handles.uitable1,'ColumnName',CLabel);
set(handles.uitable1,'Data',res);

if ~isempty(ud.t)
    [tmin,tminpos] = min(ud.t);
    [tmax,tmaxpos] = max(ud.t);
    txtout = [{'Max-Min information: '}; {sprintf('t-max: %5.3f at %s',tmax,CLabel{tmaxpos})};{sprintf('t-min: %5.3f at %s',tmin,CLabel{tminpos})};sprintf('TANOVA: %4.4f',ud.Tanova)];
    set(handles.textInfo,'String',txtout);
end

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes Ragu_TMaps wait for user response (see UIRESUME)
% uiwait(handles.figure1);


function MapBtCallback(obj, event, data,montage,style,step)

figure;
subplot('Position',[0.1 0.2 0.8 0.8]);
%for i = 1:numel(montage)
%    lbl{i} = montage(i).Name;
%end
RaguDSPMap(data,montage,style,'NoScale','Resolution',1,'Step',step);%,'Label',lbl');
subplot('Position',[0.3 0.1 0.4 0.05]);
dspCMapColorbar(step,'br');

function PMapBtCallback(obj, event, data,montage,LevelList,lbl,lblidx,llColorBar)

figure;
subplot('Position',[0.1 0.2 0.8 0.8]);
%for i = 1:numel(montage)
%    lbl{i} = montage(i).Name;
%end
RaguDSPMap(data,montage,2,'NoScale','Resolution',1,'LevelList',LevelList,'Label',lbl,'LabelIndex',lblidx);
subplot('Position',[0.3 0.1 0.4 0.05]);
dspPMapColorbar(llColorBar);





% --- Outputs from this function are returned to the command line.
function varargout = Ragu_TMaps_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function editPos_Callback(hObject, eventdata, handles)
% hObject    handle to editPos (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editPos as text
%        str2double(get(hObject,'String')) returns contents of editPos as a double


% --- Executes during object creation, after setting all properties.
function editPos_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editPos (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editNeg_Callback(hObject, eventdata, handles)
% hObject    handle to editNeg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editNeg as text
%        str2double(get(hObject,'String')) returns contents of editNeg as a double


% --- Executes during object creation, after setting all properties.
function editNeg_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editNeg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editT_Callback(hObject, eventdata, handles)
% hObject    handle to editT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editT as text
%        str2double(get(hObject,'String')) returns contents of editT as a double


% --- Executes during object creation, after setting all properties.
function editT_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
