function varargout = Ragu_ImportIndFeatFromExcel(varargin)
% RAGU_IMPORTINDFEATFROMEXCEL MATLAB code for Ragu_ImportIndFeatFromExcel.fig
%      RAGU_IMPORTINDFEATFROMEXCEL, by itself, creates a new RAGU_IMPORTINDFEATFROMEXCEL or raises the existing
%      singleton*.
%
%      H = RAGU_IMPORTINDFEATFROMEXCEL returns the handle to a new RAGU_IMPORTINDFEATFROMEXCEL or the handle to
%      the existing singleton*.
%
%      RAGU_IMPORTINDFEATFROMEXCEL('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in RAGU_IMPORTINDFEATFROMEXCEL.M with the given input arguments.
%
%      RAGU_IMPORTINDFEATFROMEXCEL('Property','Value',...) creates a new RAGU_IMPORTINDFEATFROMEXCEL or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Ragu_ImportIndFeatFromExcel_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Ragu_ImportIndFeatFromExcel_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Ragu_ImportIndFeatFromExcel

% Last Modified by GUIDE v2.5 14-Jun-2016 10:17:54

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Ragu_ImportIndFeatFromExcel_OpeningFcn, ...
                   'gui_OutputFcn',  @Ragu_ImportIndFeatFromExcel_OutputFcn, ...
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


% --- Executes just before Ragu_ImportIndFeatFromExcel is made visible.
function Ragu_ImportIndFeatFromExcel_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Ragu_ImportIndFeatFromExcel (see VARARGIN)

% Choose default command line output for Ragu_ImportIndFeatFromExcel
handles.output = hObject;

global RaguData
RaguData = varargin{4};
global nSID
nSID = str2double(get(handles.nSID,'String'));
CheckIfWeAreReady(handles);


% Update handles structure
guidata(hObject, handles);

% UIWAIT makes Ragu_ImportIndFeatFromExcel wait for user response (see UIRESUME)
uiwait(handles.Ragu_ImportIndFeatFromExcel);


% --- Outputs from this function are returned to the command line.
function varargout = Ragu_ImportIndFeatFromExcel_OutputFcn(hObject, eventdata, handles) 
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



function EditExcelFile_Callback(hObject, eventdata, handles)
% hObject    handle to EditExcelFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of EditExcelFile as text
%        str2double(get(hObject,'String')) returns contents of EditExcelFile as a double


% --- Executes during object creation, after setting all properties.
function EditExcelFile_CreateFcn(hObject, eventdata, handles)
% hObject    handle to EditExcelFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in ButtonXLSFile.
function ButtonXLSFile_Callback(hObject, eventdata, handles)
% hObject    handle to ButtonXLSFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[fnx,pnx] = uigetfile('*.xlsx','Name of the Excel data file');

if fnx == 0
    return
end

set(handles.EditExcelFile,'String',fullfile(pnx,fnx));

Sheet = inputdlg('Enter sheet name (or cancel to select the first sheet)');

if isempty(Sheet)
    Sheet = 1;
end

global xlsraw

[~,~,xlsraw] = xlsread(fullfile(pnx,fnx),Sheet);

set(handles.PopupSID,'String',xlsraw(1,:));
set(handles.PopupVID,'String',xlsraw(1,:));

global SID
SID = 1;

global VID
VID = 1;

CheckIfWeAreReady(handles);



% --- Executes on selection change in PopupSID.
function PopupSID_Callback(hObject, eventdata, handles)
% hObject    handle to PopupSID (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns PopupSID contents as cell array
%        contents{get(hObject,'Value')} returns selected item from PopupSID


global SID
SID = get(hObject,'Value');
CheckIfWeAreReady(handles);



% --- Executes during object creation, after setting all properties.
function PopupSID_CreateFcn(hObject, eventdata, handles)
% hObject    handle to PopupSID (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in PopupVID.
function PopupVID_Callback(hObject, eventdata, handles)
% hObject    handle to PopupVID (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns PopupVID contents as cell array
%        contents{get(hObject,'Value')} returns selected item from PopupVID


global VID
VID = get(hObject,'Value');
CheckIfWeAreReady(handles);

function CheckIfWeAreReady(handles)

global RaguData
global xlsraw
global SID
global VID
global nSID
WeAreReady = ~isempty(RaguData) & ~isempty(xlsraw) & ~isempty(SID) & ~isempty(VID);

if WeAreReady
    set(handles.ButtonSave,'Enable','On');
else
    set(handles.ButtonSave,'Enable','Off');
end




% --- Executes during object creation, after setting all properties.
function PopupVID_CreateFcn(hObject, eventdata, handles)
% hObject    handle to PopupVID (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in ButtonSave.
function ButtonSave_Callback(hObject, eventdata, handles)
% hObject    handle to ButtonSave (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global RaguData
global xlsraw
global SID
global VID
global nSID

rd = RaguData;
ExcludedCases = isnan(rd.IndFeature);
rd.IndFeature = nan(size(rd.IndFeature));

nWarning = 0;
Warning{1} = 'Excluded cases:';

for s = 1:size(rd.Names,1)
    idx = find(strncmpi(rd.Names(s,1),xlsraw(2:end,SID),nSID));
    switch numel(idx)
        case 1
            rd.IndFeature(s) = xlsraw{idx+1,VID};
        case 0
            nWarning = nWarning + 1;
            Warning{nWarning+1} = rd.Names{s,1};
        otherwise
            error(['Doublicate data found for ERP ' rd.Names{s,1}]);
    end
end

rd.IndFeature(ExcludedCases) = nan;
rd.IndName = xlsraw{1,VID};


if nWarning > 0
    uiwait(msgbox(Warning,'Done','modal'));
else
    uiwait(msgbox('Done','modal'));
end

switch questdlg('The chosen variable is: ','Define variable type','Continuous / rank scaled','Categorical','Categorical');
    case 'Continuous / rank scaled'
        rd.ContBetween = 1;
    case 'Categorical'
        rd.ContBetwwen = 0;
end

set(handles.output,'UserData',rd);

uiresume(handles.Ragu_ImportIndFeatFromExcel);


function nSID_Callback(hObject, eventdata, handles)
% hObject    handle to nSID (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of nSID as text
%        str2double(get(hObject,'String')) returns contents of nSID as a double

global nSID
nSID = str2double(get(hObject,'String'));

% --- Executes during object creation, after setting all properties.
function nSID_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nSID (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
