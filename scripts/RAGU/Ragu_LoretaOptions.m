function varargout = Ragu_LoretaOptions(varargin)
% RAGU_LORETAOPTIONS M-file for Ragu_LoretaOptions.fig
%      RAGU_LORETAOPTIONS, by itself, creates a new RAGU_LORETAOPTIONS or raises the existing
%      singleton*.
%
%      H = RAGU_LORETAOPTIONS returns the handle to a new RAGU_LORETAOPTIONS or the handle to
%      the existing singleton*.
%
%      RAGU_LORETAOPTIONS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in RAGU_LORETAOPTIONS.M with the given input arguments.
%
%      RAGU_LORETAOPTIONS('Property','Value',...) creates a new RAGU_LORETAOPTIONS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Ragu_LoretaOptions_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Ragu_LoretaOptions_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Ragu_LoretaOptions

% Last Modified by GUIDE v2.5 19-May-2011 07:42:53

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Ragu_LoretaOptions_OpeningFcn, ...
                   'gui_OutputFcn',  @Ragu_LoretaOptions_OutputFcn, ...
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


% --- Executes just before Ragu_LoretaOptions is made visible.
function Ragu_LoretaOptions_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Ragu_LoretaOptions (see VARARGIN)

% Choose default command line output for Ragu_LoretaOptions
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);
out = varargin{4};

if ~isfield(out,'LoretaSysDir')
    out.LoretaSysDir = '';
end

if ~isfield(out,'LoretaSPINV')
    out.LoretaSPINV = '';
end

set(handles.syspath,'String',out.LoretaSysDir);
set(handles.spinv  ,'String',out.LoretaSPINV);

% UIWAIT makes Ragu_LoretaOptions wait for user response (see UIRESUME)
uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Ragu_LoretaOptions_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
if isempty(handles)
    varargout{1} = [];
else
    varargout{1} = handles.output;
end



function syspath_Callback(hObject, eventdata, handles)
% hObject    handle to syspath (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of syspath as text
%        str2double(get(hObject,'String')) returns contents of syspath as a double


% --- Executes during object creation, after setting all properties.
function syspath_CreateFcn(hObject, eventdata, handles)
% hObject    handle to syspath (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function spinv_Callback(hObject, eventdata, handles)
% hObject    handle to spinv (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of spinv as text
%        str2double(get(hObject,'String')) returns contents of spinv as a double


% --- Executes during object creation, after setting all properties.
function spinv_CreateFcn(hObject, eventdata, handles)
% hObject    handle to spinv (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Ok.
function Ok_Callback(hObject, eventdata, handles)
% hObject    handle to Ok (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

ud.LoretaSysDir = get(handles.syspath,'String');
ud.LoretaSPINV  = get(handles.spinv  ,'String');

set(handles.output,'UserData',ud);
uiresume(handles.figure1);


% --- Executes on button press in BrowseSysPath.
function BrowseSysPath_Callback(hObject, eventdata, handles)
% hObject    handle to BrowseSysPath (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

DirName = uigetdir(get(handles.syspath,'String'),'Select directory with the ERP data');
if DirName ~= 0
    set(handles.syspath,'String',DirName);
end


% --- Executes on button press in browseSPINV.
function browseSPINV_Callback(hObject, eventdata, handles)
% hObject    handle to browseSPINV (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[filename, pathname] = uigetfile({'*.spinv','sLORETA SPINV-File';}, 'Load sLORETA SPINV File');

if isequal(filename,0) || isequal(pathname,0)
    return
else
    set(handles.spinv,'String',fullfile(pathname, filename));
end