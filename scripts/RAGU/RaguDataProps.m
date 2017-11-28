function varargout = RaguDataProps(varargin)
% RAGUDATAPROPS M-file for RaguDataProps.fig
%      RAGUDATAPROPS, by itself, creates a new RAGUDATAPROPS or raises the existing
%      singleton*.
%
%      H = RAGUDATAPROPS returns the handle to a new RAGUDATAPROPS or the handle to
%      the existing singleton*.
%
%      RAGUDATAPROPS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in RAGUDATAPROPS.M with the given input arguments.
%
%      RAGUDATAPROPS('Property','Value',...) creates a new RAGUDATAPROPS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before RaguDataProps_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to RaguDataProps_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help RaguDataProps

% Copyright 2009-2011 Thomas Koenig
% distributed under the terms of the GNU AFFERO General Public License

% Last Modified by GUIDE v2.5 07-Mar-2012 16:57:24

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @RaguDataProps_OpeningFcn, ...
                   'gui_OutputFcn',  @RaguDataProps_OutputFcn, ...
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


% --- Executes just before RaguDataProps is made visible.
function RaguDataProps_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to RaguDataProps (see VARARGIN)

% Choose default command line output for RaguDataProps
handles.output = hObject;
out = varargin{4};
set(handles.Samples_Per_Sec,'String',num2str(out.DeltaX));
set(handles.TimeOnset,'String',num2str(out.TimeOnset));

if isfield(out,'txtX')
    set(handles.txtX,'String',out.txtX);
end
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes RaguDataProps wait for user response (see UIRESUME)
uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = RaguDataProps_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure

if (~isempty(handles))
    varargout{1} = handles.output;
else
    varargout{1} = [];
end


function TimeOnset_Callback(hObject, eventdata, handles)
% hObject    handle to TimeOnset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of TimeOnset as text
%        str2double(get(hObject,'String')) returns contents of TimeOnset as a double


% --- Executes during object creation, after setting all properties.
function TimeOnset_CreateFcn(hObject, eventdata, handles)
% hObject    handle to TimeOnset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end





function Samples_Per_Sec_Callback(hObject, eventdata, handles)
% hObject    handle to Samples_Per_Sec (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Samples_Per_Sec as text
%        str2double(get(hObject,'String')) returns contents of Samples_Per_Sec as a double


% --- Executes during object creation, after setting all properties.
function Samples_Per_Sec_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Samples_Per_Sec (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in DoneButton.
function DoneButton_Callback(hObject, eventdata, handles)
% hObject    handle to DoneButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
out.TimeOnset = str2double(get(handles.TimeOnset,'String'));
out.DeltaX = str2double(get(handles.Samples_Per_Sec,'String'));
out.txtX = get(handles.txtX,'String');

set(handles.output,'UserData',out);

uiresume(handles.figure1);



function txtX_Callback(hObject, eventdata, handles)
% hObject    handle to txtX (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtX as text
%        str2double(get(hObject,'String')) returns contents of txtX as a double


% --- Executes during object creation, after setting all properties.
function txtX_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtX (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
