function varargout = Randomizer_Tanova_Options(varargin)
% RANDOMIZER_TANOVA_OPTIONS M-file for Randomizer_Tanova_Options.fig
%      RANDOMIZER_TANOVA_OPTIONS, by itself, creates a new RANDOMIZER_TANOVA_OPTIONS or raises the existing
%      singleton*.
%
%      H = RANDOMIZER_TANOVA_OPTIONS returns the handle to a new RANDOMIZER_TANOVA_OPTIONS or the handle to
%      the existing singleton*.
%
%      RANDOMIZER_TANOVA_OPTIONS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in RANDOMIZER_TANOVA_OPTIONS.M with the given input arguments.
%
%      RANDOMIZER_TANOVA_OPTIONS('Property','Value',...) creates a new RANDOMIZER_TANOVA_OPTIONS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Randomizer_Tanova_Options_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Randomizer_Tanova_Options_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Copyright 2009-2011 Thomas Koenig
% distributed under the terms of the GNU AFFERO General Public License

% Edit the above text to modify the response to help Randomizer_Tanova_Options

% Last Modified by GUIDE v2.5 04-May-2009 14:28:59

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Randomizer_Tanova_Options_OpeningFcn, ...
                   'gui_OutputFcn',  @Randomizer_Tanova_Options_OutputFcn, ...
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


% --- Executes just before Randomizer_Tanova_Options is made visible.
function Randomizer_Tanova_Options_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Randomizer_Tanova_Options (see VARARGIN)

% Choose default command line output for Randomizer_Tanova_Options
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);
out = varargin{4};
NoMean = varargin{5};

NoChangeWindow = false;
if numel(varargin) > 5
    if varargin{6} == true
        NoChangeWindow = true;
    end
end

if ~isfield(out,'StartFrame')
    out.StartFrame = [];
end

if ~isfield(out,'MeanInterval')
    out.MeanInterval = 0;
end

if ~isfield(out,'EndFrame')
    out.EndFrame = [];
end

if isempty(out.StartFrame);
    out.StartFrame = 1;
end

if isempty(out.EndFrame);
    out.EndFrame = size(out.V,4);
end

if NoMean == true
    set(handles.CBAverage,'Visible','off');
end

PreviousValues.OldStartFrame = out.StartFrame;
PreviousValues.OldEndFrame = out.EndFrame;
PreviousValues.OldMeanInterval = out.MeanInterval;
PreviousValues.NoChangeWindow = NoChangeWindow;


set(handles.EditInterValStart,'String',sprintf('%3.1f',(out.StartFrame-1)*out.DeltaX + out.TimeOnset));
set(handles.EditInterValEnd  ,'String',sprintf('%3.1f',(out.EndFrame  -1)*out.DeltaX + out.TimeOnset));
set(handles.text_IS,'String',['Interval Start (' out.txtX ')']);
set(handles.text_IE,'String',['Interval End (' out.txtX ')']);
set(handles.EditInterValStart,'UserData',PreviousValues);
set(handles.CBAverage,'Value',out.MeanInterval);
set(handles.output,'UserData',out);

% UIWAIT makes Randomizer_Tanova_Options wait for user response (see UIRESUME)
uiwait(handles.figure1);

% --- Outputs from this function are returned to the command line.
function varargout = Randomizer_Tanova_Options_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure

if (isempty(handles))
    varargout{1} = [];
else
    varargout{1} = handles.output;
end



function EditInterValStart_Callback(hObject, eventdata, handles)
% hObject    handle to EditInterValStart (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of EditInterValStart as text
%        str2double(get(hObject,'String')) returns contents of EditInterValStart as a double
out = get(handles.output,'UserData');
af  = round((str2double(get(handles.EditInterValStart,  'String'))-out.TimeOnset) /  out.DeltaX) +1;
out.StartFrame   = 1;
if af < 1
    msgbox('Chosen interval begins before the data and has been adjusted.');
    set(handles.EditInterValStart  ,'String',sprintf('%3.1f',(out.StartFrame  -1)*out.DeltaX + out.TimeOnset));
end



% --- Executes during object creation, after setting all properties.
function EditInterValStart_CreateFcn(hObject, eventdata, handles)
% hObject    handle to EditInterValStart (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function EditInterValEnd_Callback(hObject, eventdata, handles)
% hObject    handle to EditInterValEnd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of EditInterValEnd as text
%        str2double(get(hObject,'String')) returns contents of EditInterValEnd as a double

out = get(handles.output,'UserData');
ef  = round((str2double(get(handles.EditInterValEnd,  'String'))-out.TimeOnset) /  out.DeltaX) +1;
out.EndFrame   = size(out.V,4);
if ef > size(out.V,4)
    msgbox('Chosen interval ends behind the data and has been adjusted.');
    set(handles.EditInterValEnd  ,'String',sprintf('%3.1f',(out.EndFrame  -1)*out.DeltaX + out.TimeOnset));
end


% --- Executes during object creation, after setting all properties.
function EditInterValEnd_CreateFcn(hObject, eventdata, handles)
% hObject    handle to EditInterValEnd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in CBAverage.
function CBAverage_Callback(hObject, eventdata, handles)
% hObject    handle to CBAverage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of CBAverage


% --- Executes on button press in ButtonOK.
function ButtonOK_Callback(hObject, eventdata, handles)
% hObject    handle to ButtonOK (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

out = get(handles.output,'UserData');

NewStartFrame = round((str2double(get(handles.EditInterValStart,'String'))-out.TimeOnset) / out.DeltaX) +1;
NewEndFrame   = round((str2double(get(handles.EditInterValEnd,  'String'))-out.TimeOnset) / out.DeltaX) +1;
NewMeanInterval = get(handles.CBAverage,'Value');
PreviousValues = get(handles.EditInterValStart,'UserData');


if (PreviousValues.OldEndFrame ~= NewEndFrame || PreviousValues.OldStartFrame ~= NewStartFrame || PreviousValues.OldMeanInterval ~= NewMeanInterval) && PreviousValues.NoChangeWindow == false
    ButtonName = questdlg('The data window is changing. Proceeding will clear previous results.', 'Warning', 'Ok', 'Cancel','Ok');
    if strcmp(ButtonName,'Cancel')
        return
    else
        out.pTopCons = [];
        out.MeanGFP = [];
        
    end
end
    

out.StartFrame = NewStartFrame;
out.EndFrame = NewEndFrame;
out.MeanInterval = NewMeanInterval; 
out.MeanInterval = get(handles.CBAverage,'Value');
out.Continue   = 1;
set(handles.output,'UserData',out);

uiresume(handles.figure1);


% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
delete(hObject);




% --- Executes on button press in ButtonReset.
function ButtonReset_Callback(hObject, eventdata, handles)
% hObject    handle to ButtonReset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

out = get(handles.output,'UserData');
out.StartFrame = 1;
out.EndFrame   = size(out.V,4);
out.MeanInterval = 0;

set(handles.EditInterValStart,'String',sprintf('%3.1f',(out.StartFrame-1) * out.DeltaX + out.TimeOnset));
set(handles.EditInterValEnd  ,'String',sprintf('%3.1f',(out.EndFrame  -1) * out.DeltaX + out.TimeOnset));
set(handles.text_IS,'String',['Interval Start (' out.txtX ')']);
set(handles.text_IE,'String',['Interval End (' out.txtX ')']);
    
set(handles.CBAverage,'Value',out.MeanInterval);
set(handles.output,'UserData',out);

