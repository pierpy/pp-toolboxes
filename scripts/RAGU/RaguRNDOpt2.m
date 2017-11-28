function varargout = RaguRNDOpt2(varargin)
% RAGURNDOPT2 M-file for RaguRNDOpt2.fig
%      RAGURNDOPT2, by itself, creates a new RAGURNDOPT2 or raises the existing
%      singleton*.
%
%      H = RAGURNDOPT2 returns the handle to a new RAGURNDOPT2 or the handle to
%      the existing singleton*.
%
%      RAGURNDOPT2('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in RAGURNDOPT2.M with the given input arguments.
%
%      RAGURNDOPT2('Property','Value',...) creates a new RAGURNDOPT2 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before RaguRNDOpt2_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to RaguRNDOpt2_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help RaguRNDOpt2


% Copyright 2009-2011 Thomas Koenig
% distributed under the terms of the GNU AFFERO General Public License

% Last Modified by GUIDE v2.5 17-Aug-2010 13:11:52

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @RaguRNDOpt2_OpeningFcn, ...
                   'gui_OutputFcn',  @RaguRNDOpt2_OutputFcn, ...
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


% --- Executes just before RaguRNDOpt2 is made visible.
function RaguRNDOpt2_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to RaguRNDOpt2 (see VARARGIN)

% Choose default command line output for RaguRNDOpt2
handles.output = hObject;

in = varargin{4};

out.Threshold = in.Threshold;
out.NoXing = in.NoXing;
out.Iterations = in.Iterations;
out.Normalize = in.Normalize;
set(handles.output,'UserData',out);

set(handles.nIter,'String',num2str(out.Iterations,'%1.0f'));
set(handles.Threshold,'String',num2str(out.Threshold));
set(handles.NoXing,'Value',out.NoXing);
set(handles.NormListBox,'Value',out.Normalize);

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes RaguRNDOpt2 wait for user response (see UIRESUME)
uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = RaguRNDOpt2_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure

if ~isempty(handles)
    varargout{1} = handles.output;
else 
    varargout{1} = [];
end

% --- Executes on button press in DoneBT.
function DoneBT_Callback(hObject, eventdata, handles)
% hObject    handle to DoneBT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
out.Threshold = str2double(get(handles.Threshold,'String'));
out.NoXing = get(handles.NoXing,'Value');
out.Iterations = str2double(get(handles.nIter,'String'));
out.Normalize = get(handles.NormListBox,'Value');
set(handles.output,'UserData',out);
uiresume(handles.figure1);


function nIter_Callback(hObject, eventdata, handles)
% hObject    handle to nIter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of nIter as text
%        str2double(get(hObject,'String')) returns contents of nIter as a double


% --- Executes during object creation, after setting all properties.
function nIter_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nIter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Threshold_Callback(hObject, eventdata, handles)
% hObject    handle to Threshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Threshold as text
%        str2double(get(hObject,'String')) returns contents of Threshold as a double


% --- Executes during object creation, after setting all properties.
function Threshold_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Threshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in NoXing.
function NoXing_Callback(hObject, eventdata, handles)
% hObject    handle to NoXing (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of NoXing


% --- Executes on selection change in NormListBox.
function NormListBox_Callback(hObject, eventdata, handles)
% hObject    handle to NormListBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns NormListBox contents as cell array
%        contents{get(hObject,'Value')} returns selected item from NormListBox


% --- Executes during object creation, after setting all properties.
function NormListBox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to NormListBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
