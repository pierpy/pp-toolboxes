function varargout = Randomizer_Import_ERP(varargin)
% RANDOMIZER_IMPORT_ERP M-file for Randomizer_Import_ERP.fig
%      RANDOMIZER_IMPORT_ERP, by itself, creates a new RANDOMIZER_IMPORT_ERP or raises the existing
%      singleton*.
%
%      H = RANDOMIZER_IMPORT_ERP returns the handle to a new RANDOMIZER_IMPORT_ERP or the handle to
%      the existing singleton*.
%
%      RANDOMIZER_IMPORT_ERP('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in RANDOMIZER_IMPORT_ERP.M with the given input arguments.
%
%      RANDOMIZER_IMPORT_ERP('Property','Value',...) creates a new RANDOMIZER_IMPORT_ERP or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Randomizer_Import_ERP_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Randomizer_Import_ERP_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Randomizer_Import_ERP

% Copyright 2009-2011 Thomas Koenig
% distributed under the terms of the GNU AFFERO General Public License

% Last Modified by GUIDE v2.5 01-Sep-2010 11:09:34

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Randomizer_Import_ERP_OpeningFcn, ...
                   'gui_OutputFcn',  @Randomizer_Import_ERP_OutputFcn, ...
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



% --- Executes just before Randomizer_Import_ERP is made visible.
function Randomizer_Import_ERP_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Randomizer_Import_ERP (see VARARGIN)

% Choose default command line output for Randomizer_Import_ERP
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);
out = varargin{4};

%out = get(handles.output,'UserData');

if ~isfield(out,'Directory')
    out.Directory = '';
end

if ~isfield(out,'conds')
    out.conds = [];
end

if ~isfield(out,'MakeAR')
    out.MakeAR = 1;
end

if ~isfield(out,'Mask')
    out.Mask = '';
end

out.V = [];
out.Names = '';
%out.conds = [];

set(handles.DirName,'String',out.Directory);
set(handles.Mask,'String',out.Mask);
set(handles.CondList,'String',out.conds);
set(handles.CondList,'String',out.conds);
set(handles.CBAverageReference,'Value',out.MakeAR);

% UIWAIT makes Randomizer_Import_ERP wait for user response (see UIRESUME)
uiwait(handles.ImportERPDlg);

% --- Outputs from this function are returned to the command line.
function varargout = Randomizer_Import_ERP_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% Get default command line output from handles structure

if ~isempty(handles);
    varargout{1} = handles.output;
else
    varargout{1} = [];
end
% --- Executes on button press in BrowseDir.
function BrowseDir_Callback(hObject, eventdata, handles)
% hObject    handle to BrowseDir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

DirName = uigetdir(get(handles.DirName,'String'),'Select directory with the ERP data');
if DirName ~= 0
    set(handles.DirName,'String',DirName);
end

function DirName_Callback(hObject, eventdata, handles)
% hObject    handle to DirName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of DirName as text
%        str2double(get(hObject,'String')) returns contents of DirName as a double


% --- Executes during object creation, after setting all properties.
function DirName_CreateFcn(hObject, eventdata, handles)
% hObject    handle to DirName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
set(hObject,'String',pwd);

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in CondList.
function CondList_Callback(hObject, eventdata, handles)
% hObject    handle to CondList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns CondList contents as cell array
%        contents{get(hObject,'Value')} returns selected item from CondList


% --- Executes during object creation, after setting all properties.
function CondList_CreateFcn(hObject, eventdata, handles)
% hObject    handle to CondList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function NewCondition_Callback(hObject, eventdata, handles)
% hObject    handle to NewCondition (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of NewCondition as text
%        str2double(get(hObject,'String')) returns contents of NewCondition as a double


% --- Executes during object creation, after setting all properties.
function NewCondition_CreateFcn(hObject, eventdata, handles)
% hObject    handle to NewCondition (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in ButtonAddCondition.
function ButtonAddCondition_Callback(hObject, eventdata, handles)
% hObject    handle to ButtonAddCondition (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if (numel(get(handles.CondList,'String'))== 0 && get(handles.CondList,'Value') == 1)
    set(handles.CondList,'String',{get(handles.NewCondition,'String')});
else
    set(handles.CondList,'String',[get(handles.CondList,'String');{get(handles.NewCondition,'String')}]);
end

% --- Executes on button press in ButtonRemove.
function ButtonRemove_Callback(hObject, eventdata, handles)
% hObject    handle to ButtonRemove (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

str = get(handles.CondList,'String');

if (numel(str) == 0)
    return;
end

v = get(handles.CondList,'Value');
str(get(handles.CondList,'Value')) = [];
if v > numel(str)
    v = numel(str);
end

if numel(str) == 0
    str = '';
    v = 1;
end

set(handles.CondList,'String',str);
set(handles.CondList,'Value',v);

% --- Executes on button press in ButtonClear.
function ButtonClear_Callback(hObject, eventdata, handles)
% hObject    handle to ButtonClear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.CondList,'String','');
set(handles.CondList,'Value',1);

% --- Executes on button press in ButtonGo.
function ButtonGo_Callback(hObject, eventdata, handles)
% hObject    handle to ButtonGo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%cd(get(handles.DirName,'String'));

[V,Names,msg]= SlurpData(get(handles.Mask,'String'),get(handles.CondList,'String'),get(handles.DirName,'String'),get(handles.CBTranspose,'Value'),1);

uiwait(msgbox(msg,'Ragu data import'));
if ~isempty(V)
    out = get(handles.output,'UserData');
    out.MakeAR = get(handles.CBAverageReference,'Value');
    if out.MakeAR == 1
        ar = mean(V,3);
        V = V - repmat(ar,[1 1 size(V,3),1]);
    end
    
    out.V = V;
    out.Directory = get(handles.DirName,'String');
    out.Mask    = get(handles.Mask,'String');
    out.Names = Names;
    
    out.conds = get(handles.CondList,'String');
    set(handles.output,'UserData',out);
    uiresume(handles.ImportERPDlg);
end
    


function Mask_Callback(hObject, eventdata, handles)
% hObject    handle to Mask (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Mask as text
%        str2double(get(hObject,'String')) returns contents of Mask as a double


% --- Executes during object creation, after setting all properties.
function Mask_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Mask (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in ViewDir.
function ViewDir_Callback(hObject, eventdata, handles)
% hObject    handle to ViewDir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if ispc
    winopen(get(handles.DirName,'String'));
end
if ismac
    system(['open ' get(handles.DirName,'String') ' &']);
end
    

% --- Executes on button press in CBTranspose.
function CBTranspose_Callback(hObject, eventdata, handles)
% hObject    handle to CBTranspose (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of CBTranspose


% --- Executes on button press in CBAverageReference.
function CBAverageReference_Callback(hObject, eventdata, handles)
% hObject    handle to CBAverageReference (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of CBAverageReference
out.MakeAR = get(handles.CBAverageReference,'Value');