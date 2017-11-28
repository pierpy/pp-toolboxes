function varargout = Ragu_TMapper(varargin)
% RAGU_TMAPPER M-file for Ragu_TMapper.fig
%      RAGU_TMAPPER, by itself, creates a new RAGU_TMAPPER or raises the existing
%      singleton*.
%
%      H = RAGU_TMAPPER returns the handle to a new RAGU_TMAPPER or the handle to
%      the existing singleton*.
%
%      RAGU_TMAPPER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in RAGU_TMAPPER.M with the given input arguments.
%
%      RAGU_TMAPPER('Property','Value',...) creates a new RAGU_TMAPPER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Ragu_TMapper_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Ragu_TMapper_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Ragu_TMapper

% Last Modified by GUIDE v2.5 19-Jun-2015 19:19:30

% Begin initialization code - DO NOT EDIT

gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Ragu_TMapper_OpeningFcn, ...
                   'gui_OutputFcn',  @Ragu_TMapper_OutputFcn, ...
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


% --- Executes just before Ragu_TMapper is made visible.
function Ragu_TMapper_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Ragu_TMapper (see VARARGIN)

% Choose default command line output for Ragu_TMapper
handles.output = hObject;


ud.CondLabels  = varargin{4}.conds;
ud.ContBetween = varargin{4}.ContBetween;

if ud.ContBetween == 0
    ud.GroupLabels = varargin{4}.GroupLabels;
else
    ud.GroupLabels = {'All subjects'};
end

if get(handles.UseBaseline,'Value') == 0
    DoEnable = 'off';
else
    DoEnable = 'on';
end

if ~isfield(varargin{4},'LoretaSysDir')
    ud.LoretaSysDir = [];
else
    ud.LoretaSysDir = varargin{4}.LoretaSysDir;
end

if ~isfield(varargin{4},'LoretaSPINV')
    ud.LoretaSPINV = [];
else
    ud.LoretaSPINV = varargin{4}.LoretaSPINV;
end

if ~isempty(ud.LoretaSPINV)
    ud.LorInfo = Ragu_ReadSPINV(ud.LoretaSPINV);
else
    ud.LorInfo = [];
end

set(handles.GroupPosList,'String',ud.GroupLabels,'Min',0,'Max',100);
set(handles.GroupNegList,'String',ud.GroupLabels,'Min',0,'Max',100);
set(handles.CondPosList,'String',ud.CondLabels,'Min',0,'Max',100);
set(handles.CondNegList,'String',ud.CondLabels,'Min',0,'Max',100);
set(handles.CondPosBLList,'Enable',DoEnable,'String',ud.CondLabels,'Min',0,'Max',100);
set(handles.CondNegBLList,'Enable',DoEnable,'String',ud.CondLabels,'Min',0,'Max',100);

ud.GPosIdx = 1;
ud.GNegIdx = 1;
ud.CPosIdx = 1;
ud.CNegIdx = 1;
ud.BPosIdx = 1;
ud.BNegIdx = 1;

ud.Montage = varargin{4}.Channel;
ud.TimeOnset = varargin{4}.TimeOnset;
ud.DeltaX = varargin{4}.DeltaX;
ud.V = varargin{4}.V;
ud.IndFeature = varargin{4}.IndFeature;
ud.MapStyle = varargin{4}.MapStyle;
ud.txtX = varargin{4}.txtX;
ud.Normalize = varargin{4}.Normalize > 1;

if ud.Normalize
    set(handles.Normalize,'String','Normalized');
else
    set(handles.Normalize,'String','Not normalized');
end

set(handles.TTestType,'String','none');

if ~isfield(varargin{4},'StartFrame')
    ud.StartFrame = 0;
else
    ud.StartFrame = varargin{4}.StartFrame;
end

if ~isfield(varargin{4},'EndFrame')
    ud.EndFrame = size(ud.V,4);
else
    ud.EndFrame = varargin{4}.EndFrame;
end

ud.time = (0:(size(ud.V,4)-1));
ud.time = ud.time * ud.DeltaX + ud.TimeOnset;
set(handles.text_IStart,'String',['Start ' ud.txtX]);
set(handles.text_IEnd  ,'String',['End ' ud.txtX]);

set(handles.EditStart,'String',sprintf('%3.1f',ud.time(ud.StartFrame)));
set(handles.EditEnd  ,'String',sprintf('%3.1f',ud.time(ud.EndFrame)));
    
set(handles.figure1,'UserData',ud);

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes Ragu_TMapper wait for user response (see UIRESUME)
uiwait(handles.figure1);

function TMapCellSelection(hObject, eventdata, handles,tableID)

ud = get(handles.figure1,'UserData');

%set(handles.sLoreta,'Enable','off');

switch tableID
    case 1
        ud.CPosIdx = eventdata.Indices(:,1);
    case 2  
        ud.CNegIdx = eventdata.Indices(:,1);
    case 3  
        ud.GPosIdx = eventdata.Indices(:,1);
    case 4  
        ud.GNegIdx = eventdata.Indices(:,1);
    case 5
        ud.BPosIdx = eventdata.Indices(:,1);
    case 6  
        ud.BNegIdx = eventdata.Indices(:,1);

end

set(handles.figure1,'UserData',ud);


function res = GetSubjectSelection(ud,idx)
res = zeros(size(ud.V,1),1);

if ud.ContBetween == 0
    for i = 1:numel(idx)
        res(ud.IndFeature == idx(i)) = 1;
    end
else
    res = ~isnan(ud.IndFeature);
end

function res = PreprocessedMean(d,nrm,logFlag)
%size(d)
if nargin < 3
    logFlag = 0;
end

if nrm == 1
    d = NormDimL2(d,3);
end

if nrm == 2
    m = mean(d,3);
    d = d ./ repmat(m,[1,1,size(d,3),1]);
end

if logFlag == 1
    d = log(d);
end


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

function DoTheTShow(handles,ShowTable,ShowPValues,RandRuns)

if nargin < 4
    RandRuns = 5000;
end

if nargin < 2
    ShowTable = false;
end

ud = get(handles.figure1,'UserData');
if get(handles.UseBaseline,'Value') == 1
    if isempty(ud.BPosIdx) || isempty(ud.BNegIdx)
        return
    elseif (IsIndexDifferent(ud.CPosIdx,ud.BPosIdx) == false) || (IsIndexDifferent(ud.CNegIdx,ud.BNegIdx) == false)
        return
    end
end

if isempty(ud.CPosIdx) + isempty(ud.GPosIdx) + isempty(ud.CNegIdx) + isempty(ud.GNegIdx) > 0
    return
end

TRange = str2double(get(handles.TRange,'String'));
MRange = str2double(get(handles.MeanRange,'String'));

S1ToUse = GetSubjectSelection(ud,ud.GPosIdx);
S2ToUse = GetSubjectSelection(ud,ud.GNegIdx);
    
Pos = PreprocessedMean(ud.V(S1ToUse == 1,ud.CPosIdx,:,ud.StartFrame:ud.EndFrame),ud.Normalize);
Neg = PreprocessedMean(ud.V(S2ToUse == 1,ud.CNegIdx,:,ud.StartFrame:ud.EndFrame),ud.Normalize);

if get(handles.UseBaseline,'Value') == 1
    Pos = Pos - PreprocessedMean(ud.V(S1ToUse == 1,ud.BPosIdx,:,ud.StartFrame:ud.EndFrame),ud.Normalize);
    Neg = Neg - PreprocessedMean(ud.V(S2ToUse == 1,ud.BNegIdx,:,ud.StartFrame:ud.EndFrame),ud.Normalize);
end

MPos = mean(Pos);
MNeg = mean(Neg);
        
outData.MPos = MPos;
outData.MNeg = MNeg;
outData.Montage = ud.Montage;
outData.MapStyle = ud.MapStyle;
outData.MRange   = MRange;
outData.TRange   = TRange;
outData.t  = [];
outData.ShowPValues = ShowPValues;

TTestOK = true;
if numel(ud.CPosIdx) == numel(ud.CNegIdx)
    if(sum(ud.CPosIdx ~= ud.CNegIdx) == 0)
        TTestOK = false;
    end
end

outData.PosString{1} ={'Groups:'};
idx = 1;
for i = 1:numel(ud.GPosIdx)
    outData.PosString{idx+i} = ud.GroupLabels(i);
    idx= idx+1;
end


if get(handles.UseBaseline,'Value') == 1
    outData.PosString = [{'Groups: '} ;ud.GroupLabels(ud.GPosIdx); {'Conditions: '}; ud.CondLabels(ud.CPosIdx); {'Baseline'}; ud.CondLabels(ud.BPosIdx)];
    outData.NegString = [{'Groups: '} ;ud.GroupLabels(ud.GNegIdx); {'Conditions: '}; ud.CondLabels(ud.CNegIdx); {'Baseline'}; ud.CondLabels(ud.BNegIdx)];
else
    outData.PosString = [{'Groups: '} ;ud.GroupLabels(ud.GPosIdx); {'Conditions: '}; ud.CondLabels(ud.CPosIdx)];
    outData.NegString = [{'Groups: '} ;ud.GroupLabels(ud.GNegIdx); {'Conditions: '}; ud.CondLabels(ud.CNegIdx)];
end

outData.Title = sprintf('%s (%2.1f - %2.1f)',get(handles.editTitle,'String'),ud.time(ud.StartFrame),ud.time(ud.EndFrame));

set(handles.sLoreta,'Enable','off');

if (sum(S1ToUse ~= S2ToUse) == 0)
    set(handles.TTestType,'String','paired');
    if TTestOK
        outData.t = tvalue(Pos-Neg,[],1);
        outData.Tanova = TANOVAP(Pos-Neg,RandRuns);
        outData.df = size(Pos,1) - 1;
        outData.tPos = tvalue(Pos,[],1);
%        res = [res;outData.t];
        set(handles.sLoreta,'Enable','on');
    end
else
    set(handles.TTestType,'String','unpaired');
    outData.t = tvalue(Pos,Neg,1);
    outData.Tanova = TANOVAUP(Pos,Neg,RandRuns);
    outData.df = size(Pos,1) + size(Neg,1)- 2;
    set(handles.sLoreta,'Enable','on');
end

outData.tPos = tvalue(Pos,[],1);
outData.tNeg = tvalue(Neg,[],1);

if ShowTable == false
    Ragu_TMaps([],[],[],outData);
else
    Ragu_TTable([],[],[],outData);
end

function p = TANOVAP(in,runs)
es = nan(runs,1);
es(1) = sqrt(sum(mean(in,1).^2,2));
for r = 2:runs
    perm = repmat(sign(randn(size(in,1),1)),1,size(in,2));
    es(r) = sqrt(sum(mean(in.*perm,1).^2,2));
end

p = sum(es >= es(1)) / runs;



function p = TANOVAUP(g1,g2,runs)
es = nan(runs,1);

es(1) = sqrt(sum(((mean(g1,1)-mean(g2,1)).^2),2));
%es(1) = std(mean(g1,1)-mean(g2,1));
g = [g1;g2];
for r = 2:runs
    ga = randn(size(g,1),1) > 0;
    es(r) = sqrt(sum(((mean(g(ga == true,:),1)-mean(g(ga == false,:),1)).^2),2));
end
p = sum(es >= es(1)) / runs;


% --- Outputs from this function are returned to the command line.
function varargout = Ragu_TMapper_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
if isempty(handles)
    varargout{1} = [];
else
    ud = get(handles.figure1,'UserData');
    set(handles.output,'UserData',ud);
    varargout{1} = handles.output;
end

function EditStart_Callback(hObject, eventdata, handles)
% hObject    handle to EditStart (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

out = get(handles.output,'UserData');
out.StartFrame  = round((str2double(get(handles.EditStart,  'String'))-out.TimeOnset) /  out.DeltaX) +1;

if out.StartFrame < 1
    msgbox('Chosen interval begins before the data and has been adjusted.');
    out.StartFrame = 1;
    set(handles.EditStart,'String',sprintf('%3.1f',(out.time(out.StartFrame))));
end

set(handles.figure1,'UserData',out);
TMapCellSelection([], [], handles,0);

function EditStart_CreateFcn(hObject, eventdata, handles)
% hObject    handle to EditStart (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function EditEnd_Callback(hObject, eventdata, handles)
% hObject    handle to EditEnd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of EditEnd as text
%        str2double(get(hObject,'String')) returns contents of EditEnd as a double

out = get(handles.output,'UserData');
out.EndFrame  = round((str2double(get(handles.EditEnd,  'String'))-out.TimeOnset)/ out.DeltaX) +1;

if (out.EndFrame < out.StartFrame || out.EndFrame > size(out.V,4))
    msgbox('Chosen interval end has been adjusted.');
    out.EndFrame = size(out.V,4);
    set(handles.EditEnd,'String',sprintf('%3.1f',(out.time(out.EndFrame))));
end

set(handles.figure1,'UserData',out);
TMapCellSelection([], [], handles,0);

function EditEnd_CreateFcn(hObject, eventdata, handles)
% hObject    handle to EditEnd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function Normalize_Callback(hObject, eventdata, handles)

function Normalize_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Normalize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function TTestType_Callback(hObject, eventdata, handles)

function TTestType_CreateFcn(hObject, eventdata, handles)
% hObject    handle to TTestType (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function MeanRange_Callback(hObject, eventdata, handles)
% hObject    handle to MeanRange (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%DoTheTShow(handles);

function MeanRange_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MeanRange (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function TRange_Callback(hObject, eventdata, handles)
% hObject    handle to TRange (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%DoTheTShow(handles);

function TRange_CreateFcn(hObject, eventdata, handles)
% hObject    handle to TRange (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
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

function res = V2Lor(in,lor)
s = size(in);
if numel(s) < 4
    s(4) = 1;
end
 
o = [3 1 2 4];

v = reshape(permute(in,o),[s(3),s(1)*s(2)*s(4)]);
res = Ragu_ComputeLoreta(v,lor);
res = ipermute(reshape(res,[size(lor.spinv,3),s(1),s(2),s(4)]),o);


function UseBaseline_Callback(hObject, eventdata, handles)
% hObject    handle to UseBaseline (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if get(hObject,'Value') == 0
    DoEnable = 'off';
else
    DoEnable = 'on';
end

set(handles.CondPosBLList,'Enable',DoEnable);
set(handles.CondNegBLList,'Enable',DoEnable);

TMapCellSelection(hObject, eventdata, handles,0);

function sLoreta_Callback(hObject, eventdata, handles)
% hObject    handle to sLoreta (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

ud = get(handles.figure1,'UserData');

if isempty(ud.LorInfo)
    h = Ragu_LoretaOptions(hObject, eventdata, handles,ud);
    if (isempty(h))
        return
    end
    out = get(h,'UserData');

    ud.LoretaSysDir = out.LoretaSysDir;
    ud.LoretaSPINV  = out.LoretaSPINV;
    ud.LorInfo = Ragu_ReadSPINV(ud.LoretaSPINV);
    set(handles.figure1,'UserData',ud);
    close(h);
end

S1ToUse = GetSubjectSelection(ud,ud.GPosIdx);
S2ToUse = GetSubjectSelection(ud,ud.GNegIdx);

if (sum(S1ToUse ~= S2ToUse) == 0)
    set(handles.TTestType,'String','paired');
else
    set(handles.TTestType,'String','unpaired');
end


LogFlag = 0;
if get(handles.LorLog,'Value') == 1
    LogFlag = 1;
end

s1idx = find(S1ToUse == 1);
s2idx = find(S2ToUse == 1);

Pos = zeros(numel(s1idx),ud.LorInfo.nvox);
Neg = zeros(numel(s2idx),ud.LorInfo.nvox);

for s = 1:numel(s1idx)
    if get(handles.UseBaseline,'Value') == 1
        Pos(s,:) = PreprocessedMean(V2Lor(ud.V(s1idx(s),ud.CPosIdx,:,ud.StartFrame:ud.EndFrame),ud.LorInfo),ud.Normalize*2,LogFlag) ...
                 - PreprocessedMean(V2Lor(ud.V(s1idx(s),ud.BPosIdx,:,ud.StartFrame:ud.EndFrame),ud.LorInfo),ud.Normalize*2,LogFlag);
    else
        Pos(s,:) = PreprocessedMean(V2Lor(ud.V(s1idx(s),ud.CPosIdx,:,ud.StartFrame:ud.EndFrame),ud.LorInfo),ud.Normalize*2,LogFlag);
    end
end

for s = 1:numel(s2idx)
    if get(handles.UseBaseline,'Value') == 1
        Neg(s,:) = PreprocessedMean(V2Lor(ud.V(s2idx(s),ud.CNegIdx,:,ud.StartFrame:ud.EndFrame),ud.LorInfo),ud.Normalize*2,LogFlag) ...
                 - PreprocessedMean(V2Lor(ud.V(s2idx(s),ud.BNegIdx,:,ud.StartFrame:ud.EndFrame),ud.LorInfo),ud.Normalize*2,LogFlag);
    else
        Neg(s,:) = PreprocessedMean(V2Lor(ud.V(s2idx(s),ud.CNegIdx,:,ud.StartFrame:ud.EndFrame),ud.LorInfo),ud.Normalize*2,LogFlag);    
    end
end

if (sum(S1ToUse ~= S2ToUse) == 0)
    t = tvalue(Pos-Neg,[],1);
else
    t = tvalue(Pos,Neg,1,'e');
end
li.LoretaSysDir = ud.LoretaSysDir;
li.LorImage = t;
li.Title = sprintf('%s (%2.1f - %2.1f)',get(handles.editTitle,'String'),ud.time(ud.StartFrame),ud.time(ud.EndFrame));
dspLoreta(hObject,eventdata,handles,li);

% --- Executes on button press in ButtonShowTMaps.
function ButtonShowTMaps_Callback(hObject, eventdata, handles)
% hObject    handle to ButtonShowTMaps (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

DoTheTShow(handles,false,false);

% --- Executes on button press in LorLog.
function LorLog_Callback(hObject, eventdata, handles)

% --- Executes on button press in ButtonLorOptions.
function ButtonLorOptions_Callback(hObject, eventdata, handles)
% hObject    handle to ButtonLorOptions (see GCBO)
% eventdata  reserved - to be defined in a future version of
% MATLAB

ud = get(handles.figure1,'UserData');
h = Ragu_LoretaOptions(hObject, eventdata, handles,ud);

if (isempty(h))
    return
end

out = get(h,'UserData');

ud.LoretaSysDir = out.LoretaSysDir;
ud.LoretaSPINV  = out.LoretaSPINV;
ud.LorInfo = Ragu_ReadSPINV(ud.LoretaSPINV);
close(h);

set(handles.figure1,'UserData',ud);

% --- Executes on button press in ButtonDone.
function ButtonDone_Callback(hObject, eventdata, handles)

uiresume(handles.figure1);

% --- Executes on selection change in GroupPosList.
function GroupPosList_Callback(hObject, eventdata, handles)

evt.Indices = get(hObject,'Value')';
TMapCellSelection(hObject, evt, handles,3);

% --- Executes during object creation, after setting all properties.
function GroupPosList_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in GroupNegList.
function GroupNegList_Callback(hObject, eventdata, handles)

evt.Indices = get(hObject,'Value')';
TMapCellSelection(hObject, evt, handles,4);

% --- Executes during object creation, after setting all properties.
function GroupNegList_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in CondPosList.
function CondPosList_Callback(hObject, eventdata, handles)

evt.Indices = get(hObject,'Value')';
TMapCellSelection(hObject, evt, handles,1);

% --- Executes during object creation, after setting all properties.
function CondPosList_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in CondNegList.
function CondNegList_Callback(hObject, eventdata, handles)


evt.Indices = get(hObject,'Value')';
TMapCellSelection(hObject, evt, handles,2);

% --- Executes during object creation, after setting all properties.
function CondNegList_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in CondPosBLList.
function CondPosBLList_Callback(hObject, eventdata, handles)

evt.Indices = get(hObject,'Value')';
TMapCellSelection(hObject, evt, handles,5);

% --- Executes during object creation, after setting all properties.
function CondPosBLList_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in CondNegBLList.
function CondNegBLList_Callback(hObject, eventdata, handles)
evt.Indices = get(hObject,'Value')';
TMapCellSelection(hObject, evt, handles,6);

% --- Executes during object creation, after setting all properties.
function CondNegBLList_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in ButtonShowTable.
function ButtonShowTable_Callback(hObject, eventdata, handles)

DoTheTShow(handles,true,false);

function editTitle_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function editTitle_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in ButtonShowPMaps.
function ButtonShowPMaps_Callback(hObject, eventdata, handles)
% hObject    handle to ButtonShowPMaps (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

DoTheTShow(handles,false,true);


% --- Executes on button press in checkboxTanova.
function checkboxTanova_Callback(hObject, eventdata, handles)
% hObject    handle to checkboxTanova (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkboxTanova
