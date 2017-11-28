function varargout = RaguMSParameters(varargin)
% RAGUMSPARAMETERS M-file for RaguMSParameters.fig
%      RAGUMSPARAMETERS, by itself, creates a new RAGUMSPARAMETERS or raises the existing
%      singleton*.
%
%      H = RAGUMSPARAMETERS returns the handle to a new RAGUMSPARAMETERS or the handle to
%      the existing singleton*.
%
%      RAGUMSPARAMETERS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in RAGUMSPARAMETERS.M with the given input arguments.
%
%      RAGUMSPARAMETERS('Property','Value',...) creates a new RAGUMSPARAMETERS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before RaguMSParameters_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to RaguMSParameters_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help RaguMSParameters

% Last Modified by GUIDE v2.5 14-Jun-2016 17:19:26

% Copyright 2009-2011 Thomas Koenig
% distributed under the terms of the GNU AFFERO General Public License

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @RaguMSParameters_OpeningFcn, ...
                   'gui_OutputFcn',  @RaguMSParameters_OutputFcn, ...
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


% --- Executes just before RaguMSParameters is made visible.
function RaguMSParameters_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to RaguMSParameters (see VARARGIN)

% Choose default command line output for RaguMSParameters
handles.output = hObject;

out = varargin{4};
if numel(varargin) > 4
    if varargin{5} == 1
        set(handles.nMicrostates,'Enable','off');
        set(handles.nReinitializations,'Enable','off');
        set(handles.OptStart,'enable','off');
        set(handles.OptEnd,'enable','off');
        set(handles.OptNTraining,'enable','off');
    end
end

set(handles.OptNCases,'String',num2str(sum(~isnan(out.IndFeature))));
if isfield(out,'UseAAHC')
    set(handles.checkboxAAHC,'Value',out.UseAAHC == true)
end
if isfield(out,'FixedK')
    set(handles.FixedN,'Value',out.FixedK == true);
    set(handles.XValN,'Value', out.FixedK == false);
end
if out.FixedK == true
    dummy.NewValue  = handles.FixedN;
else
    dummy.NewValue  = handles.XValN;
end

PanelOptMethod_SelectionChangeFcn(handles.PanelOptMethod, dummy, handles);

if isfield(out,'XValRestarts')
    set(handles.XValRestarts,'String',num2str(out.XValRestarts));
end

if isfield(out,'OptStart')
    set(handles.OptStart,'String',num2str(out.OptStart));
end

if isfield(out,'OptEnd')
    set(handles.OptEnd,'String',num2str(out.OptEnd));
end

if ~isfield(out,'OptNTraining')
    out.OptNTraining = floor(sum(~isnan(out.IndFeature))/2);
end
set(handles.OptNTraining,'String',num2str(out.OptNTraining));

if isfield(out,'nStates')
    set(handles.nMicrostates,'String',num2str(out.nStates));
end

if isfield(out,'nReiter')
    set(handles.nReinitializations,'String',num2str(out.nReiter));
end

if isfield(out,'bSmoothLabels')
    set(handles.cbSmoothLabels,'Value',out.bSmoothLabels);
end

if isfield(out,'nWindowSize')
    set(handles.editWindowSize,'String',num2str(out.nWindowSize));
end

if isfield(out,'LabelPenalty')
    set(handles.NSPedit,'String',num2str(out.LabelPenalty));
end

if isfield(out,'NoInconsistentMaps')
    set(handles.cbNoIncMaps,'value',out.NoInconsistentMaps);
end

if ~isfield(out,'UseAAHC')
    out.UseAAHC = false;
end
    
if out.UseAAHC
    set(handles.nReinitializations,'Enable','off');
else
    set(handles.nReinitializations,'Enable','on');
end


UpdateGUI(handles);
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes RaguMSParameters wait for user response (see UIRESUME)
uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = RaguMSParameters_OutputFcn(hObject, eventdata, handles) 
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


function nMicrostates_Callback(hObject, eventdata, handles)
% hObject    handle to nMicrostates (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of nMicrostates as text
%        str2double(get(hObject,'String')) returns contents of nMicrostates as a double


% --- Executes during object creation, after setting all properties.
function nMicrostates_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nMicrostates (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function nReinitializations_Callback(hObject, eventdata, handles)
% hObject    handle to nReinitializations (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of nReinitializations as text
%        str2double(get(hObject,'String')) returns contents of nReinitializations as a double


% --- Executes during object creation, after setting all properties.
function nReinitializations_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nReinitializations (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in ButtonOK.
function ButtonOK_Callback(hObject, eventdata, handles)
% hObject    handle to ButtonOK (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

out.nStates = round(str2double(get(handles.nMicrostates,'String')));
out.nReiter = round(str2double(get(handles.nReinitializations,'String')));
out.bSmoothLabels = get(handles.cbSmoothLabels,'Value');
out.nWindowSize = round(str2double(get(handles.editWindowSize,'String')));
out.LabelPenalty = str2double(get(handles.NSPedit,'String'));
out.NoInconsistentMaps = get(handles.cbNoIncMaps,'Value');
out.UseAAHC = get(handles.checkboxAAHC,'Value');

out.FixedK   = get(handles.FixedN,'Value');
out.OptStart = round(str2double(get(handles.OptStart,'String')));
out.OptEnd = round(str2double(get(handles.OptEnd,'String')));
out.OptNTraining = round(str2double(get(handles.OptNTraining,'String')));
out.XValRestarts = round(str2double(get(handles.XValRestarts,'String')));
set(handles.figure1,'UserData',out);
uiresume(handles.figure1);


function UpdateGUI(handles)
if get(handles.cbSmoothLabels,'Value') == 1
    enable = 'on';
else
    enable = 'off';
end

set(handles.editWindowSize,'enable',enable);
set(handles.NSPedit,'enable',enable);



% --- Executes on button press in cbSmoothLabels.
function cbSmoothLabels_Callback(hObject, eventdata, handles)
% hObject    handle to cbSmoothLabels (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
UpdateGUI(handles);

% Hint: get(hObject,'Value') returns toggle state of cbSmoothLabels



function editWindowSize_Callback(hObject, eventdata, handles)
% hObject    handle to editWindowSize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editWindowSize as text
%        str2double(get(hObject,'String')) returns contents of editWindowSize as a double


% --- Executes during object creation, after setting all properties.
function editWindowSize_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editWindowSize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function NSPedit_Callback(hObject, eventdata, handles)
% hObject    handle to NSPedit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of NSPedit as text
%        str2double(get(hObject,'String')) returns contents of NSPedit as a double


% --- Executes during object creation, after setting all properties.
function NSPedit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to NSPedit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in cbNoIncMaps.
function cbNoIncMaps_Callback(hObject, eventdata, handles)
% hObject    handle to cbNoIncMaps (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of cbNoIncMaps



function OptStart_Callback(hObject, eventdata, handles)
% hObject    handle to OptStart (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of OptStart as text
%        str2double(get(hObject,'String')) returns contents of OptStart as a double


% --- Executes during object creation, after setting all properties.
function OptStart_CreateFcn(hObject, eventdata, handles)
% hObject    handle to OptStart (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function OptEnd_Callback(hObject, eventdata, handles)
% hObject    handle to OptEnd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of OptEnd as text
%        str2double(get(hObject,'String')) returns contents of OptEnd as a double


% --- Executes during object creation, after setting all properties.
function OptEnd_CreateFcn(hObject, eventdata, handles)
% hObject    handle to OptEnd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function OptNTraining_Callback(hObject, eventdata, handles)
% hObject    handle to OptNTraining (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of OptNTraining as text
%        str2double(get(hObject,'String')) returns contents of OptNTraining as a double


% --- Executes during object creation, after setting all properties.
function OptNTraining_CreateFcn(hObject, eventdata, handles)
% hObject    handle to OptNTraining (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function OptNCases_Callback(hObject, eventdata, handles)
% hObject    handle to OptNCases (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of OptNCases as text
%        str2double(get(hObject,'String')) returns contents of OptNCases as a double


% --- Executes during object creation, after setting all properties.
function OptNCases_CreateFcn(hObject, eventdata, handles)
% hObject    handle to OptNCases (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes when selected object is changed in PanelOptMethod.
function PanelOptMethod_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in PanelOptMethod 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)

if(eventdata.NewValue == handles.FixedN)
    e1 = 'on';
    e2 = 'off';
else
    e1 = 'off';
    e2 = 'on';
end
   
set(handles.nMicrostates,'enable',e1);
set(handles.OptStart,'enable',e2);
set(handles.OptEnd,'enable',e2);
set(handles.OptNTraining,'enable',e2);



function XValRestarts_Callback(hObject, eventdata, handles)
% hObject    handle to XValRestarts (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of XValRestarts as text
%        str2double(get(hObject,'String')) returns contents of XValRestarts as a double


% --- Executes during object creation, after setting all properties.
function XValRestarts_CreateFcn(hObject, eventdata, handles)
% hObject    handle to XValRestarts (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkboxAAHC.
function checkboxAAHC_Callback(hObject, eventdata, handles)
% hObject    handle to checkboxAAHC (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkboxAAHC
if get(hObject,'Value') > 0
    set(handles.nReinitializations,'Enable','off');
else
    set(handles.nReinitializations,'Enable','on');
end
