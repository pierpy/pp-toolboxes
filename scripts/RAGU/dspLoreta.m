function varargout = dspLoreta(varargin)
% DSPLORETA M-file for dspLoreta.fig
%      DSPLORETA, by itself, creates a new DSPLORETA or raises the existing
%      singleton*.
%
%      H = DSPLORETA returns the handle to a new DSPLORETA or the handle to
%      the existing singleton*.
%
%      DSPLORETA('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in DSPLORETA.M with the given input arguments.
%
%      DSPLORETA('Property','Value',...) creates a new DSPLORETA or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before dspLoreta_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to dspLoreta_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help dspLoreta

% Last Modified by GUIDE v2.5 20-Nov-2013 16:56:48

% Begin initialization code - DO NOT EDIT
gui_Singleton = 0;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @dspLoreta_OpeningFcn, ...
                   'gui_OutputFcn',  @dspLoreta_OutputFcn, ...
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


% --- Executes just before dspLoreta is made visible.
function dspLoreta_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to dspLoreta (see VARARGIN)

% Choose default command line output for dspLoreta
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

in = varargin{4};

rawDir = dir([in.LoretaSysDir '\*.raw']);
if isempty(rawDir)
    errordlg('Loreta system directory incorrect.', 'dspLoreta');
    close(hObject);
    return;
end


if (isfield(in,'Title'))
    set(hObject,'Name',in.Title);
end

AnatomyFile = cell(numel(rawDir),1);
for i = 1:numel(rawDir)
    AnatomyFile{i} = rawDir(i).name;
end

DefaultAnatomy = 3;

set(handles.lbAnatomies,'String',AnatomyFile,'Value',DefaultAnatomy);
set(handles.lbAnatomies,'UserData',in.LoretaSysDir);

fp = fopen([in.LoretaSysDir '\' AnatomyFile{DefaultAnatomy}]);
ud.Anatomy = reshape(fread(fp,'uint8'),[181,217,181]) / 2;
fclose(fp);

fid = fopen([in.LoretaSysDir '\MNI-BAs-6239-voxels.csv'],'rt');
c = textscan(fid,'%f%f%f%s%s%s','Delimiter',',');
fclose(fid);
ud.VoxInfo = c;

ud.XOff =  91;
ud.YOff = 130;
ud.ZOff =  72;

ud.VoxLR = c{1} + ud.XOff;
ud.VoxAP = c{2} + ud.YOff;
ud.VoxSI = c{3} + ud.ZOff;
ud.Lobe = c{4};
ud.Gyrus = c{5};
ud.BA = c{6};
ud.Lori = in.LorImage;

ud.LRPos = ud.XOff;
ud.APPos = ud.YOff;
ud.SIPos = ud.ZOff;

UpdatePositionInfo(ud,handles);
ShowTalInfo(ud,handles);

ud.SVoxInd = sub2ind([29 34 24],c{1}/5 + 15,c{2}/5 +21,c{3}/5+10);
ud.V = nan(29,34,24);
ud.V(ud.SVoxInd) = in.LorImage;

[ud.SVoxGridAP,ud.SVoxGridLR,ud.SVoxGridSI] = meshgrid((1:5:166)+29,(1:5:141)+20,(1:5:116)+26);

%LR * AP * SI
ud.Threshold = str2double(get(handles.editThreshold,'String'));
ud.Max       = str2double(get(handles.editMax      ,'String'));

ud.TopViewHandle = ShowTopView(ud,handles);
hold on
ud.TopViewHC = plot([1 181],[ud.APPos ud.APPos],'-g','ButtonDownFcn',{@TopViewKlicked,handles});
ud.TopViewVC = plot([ud.LRPos ud.LRPos],[1 217],'-g','ButtonDownFcn',{@TopViewKlicked,handles});

ud.BackViewHandle = ShowBackView(ud,handles);
hold on
ud.BackViewHC = plot([1 181],[ud.SIPos ud.SIPos],'-g','ButtonDownFcn',{@BackViewKlicked,handles});
ud.BackViewVC = plot([ud.LRPos ud.LRPos],[1 181],'-g','ButtonDownFcn',{@BackViewKlicked,handles});

ud.SideViewHandle = ShowSideView(ud,handles);
hold on
ud.SideViewHC = plot([1 217],[ud.SIPos ud.SIPos],'-g','ButtonDownFcn',{@SideViewKlicked,handles});
ud.SideViewVC = plot([ud.APPos ud.APPos],[1 181],'-g','ButtonDownFcn',{@SideViewKlicked,handles});

set(handles.SideView,'YTick',22:25:172,'YTickLabel',{'-50','-25','0','25','50','75','100','100'});
set(handles.BackView,'YTick',22:25:172,'YTickLabel',{'-50','-25','0','25','50','75','100','100'});

set(handles.BackView,'XTick',16:25:166,'XTickLabel',{'-75','-50','-25','0','25','50','75'});
set(handles.TopView ,'XTick',16:25:166,'XTickLabel',{'-75','-50','-25','0','25','50','75'});

set(handles.SideView,'XTick', 5:25:205,'XTickLabel',{'-125','-100','-75','-50','-25','0','25','50','75'});
set(handles.TopView ,'YTick', 5:25:205,'YTickLabel',{'-125','-100','-75','-50','-25','0','25','50','75'});



%colormap([gray(128);BlueRed(128)]);
colormap([(1-gray(128));LorBlueRed(128)]);
ud.colormap = [(1-gray(128));LorBlueRed(128)];

ShowColorScale(handles,ud);

set(handles.figure1,'UserData',ud);

% UIWAIT makes dspLoreta wait for user response (see UIRESUME)
%uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = dspLoreta_OutputFcn(hObject, eventdata, handles) 
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

function h = ShowSideView(ud,handles)
%--------------------------------
subplot(handles.SideView);
h = imagesc(GetSideImage(ud,ud.LRPos),[0 255]);
set(h,'ButtonDownFcn',{@SideViewKlicked,handles});
axis xy

function h = ShowBackView(ud,handles)
%--------------------------------
subplot(handles.BackView);
h = imagesc(GetBackImage(ud,ud.APPos),[0 255]);
set(h,'ButtonDownFcn',{@BackViewKlicked,handles});
axis xy


function h = ShowTopView(ud,handles)
%-------------------------------
subplot(handles.TopView);
h = imagesc(GetTopImage(ud,ud.SIPos),[0 255]);
set(h,'ButtonDownFcn',{@TopViewKlicked,handles});
axis xy

function out = MyClip(in,crit)

out = in;
out(out >  crit) =  crit;
out(out < -crit) = -crit;


function res = GetTopImage(ud,pos)
%---------------------------------
%size(ud.SVoxGridAP)
%size(ud.SVoxGridLR)
%size(ud.SVoxGridSI)
%size(ud.V)

[X,Y,Z] = meshgrid(1:217,1:181,pos);

%v = interp3(ud.SVoxGridAP,ud.SVoxGridLR,ud.SVoxGridSI,ud.V,1:217,1:181,pos,'nearest')';
v = interp3(ud.SVoxGridAP,ud.SVoxGridLR,ud.SVoxGridSI,ud.V,X,Y,Z,'nearest')';
idx = find(abs(v) > ud.Threshold);
%idx(rem(idx,2) == 0) = [];
v = v / ud.Max * 63;
v = MyClip(v,63);
res = squeeze(ud.Anatomy(:,:,pos))';
res(idx) = v(idx)+196;

function res = GetBackImage(ud,pos)
%----------------------------------
v = squeeze(interp3(ud.SVoxGridAP,ud.SVoxGridLR,ud.SVoxGridSI,ud.V,pos,1:181,1:181,'nearest'))';
idx = find(abs(v) > ud.Threshold);
%idx(rem(idx,2) == 0) = [];
v = v / ud.Max * 63;
v = MyClip(v,63);
res = squeeze(ud.Anatomy(:,pos,:))';
res(idx) = v(idx)+196;

function res = GetSideImage(ud,pos)
%---------------------------------
v = squeeze(interp3(ud.SVoxGridAP,ud.SVoxGridLR,ud.SVoxGridSI,ud.V,1:217,pos,1:181,'nearest'))';
idx = find(abs(v) > ud.Threshold);
%idx(rem(idx,2) == 0) = [];
v = v / ud.Max * 63;
v = MyClip(v,63);
res = squeeze(ud.Anatomy(pos,:,:))';
res(idx) = v(idx)+196;



function SideViewKlicked(obj,event,handles)
%----------------------------------
pt = get(get(obj,'Parent'),'CurrentPoint');
ud = get(gcf,'UserData');
ud.SIPos = round(pt(1,2));
ud.APPos = round(pt(1,1));

set(ud.SideViewHC,'YData',[ud.SIPos ud.SIPos]);
set(ud.SideViewVC,'XData',[ud.APPos ud.APPos]);

set(ud.TopViewHandle,'CData',GetTopImage(ud,ud.SIPos));
set(ud.TopViewHC,'YData',[ud.APPos ud.APPos]);

set(ud.BackViewHandle,'CData',GetBackImage(ud,ud.APPos));
set(ud.BackViewHC,'YData',[ud.SIPos ud.SIPos]);

UpdatePositionInfo(ud,handles);
ShowTalInfo(ud,handles);
set(gcf,'UserData',ud);

function BackViewKlicked(obj,event,handles)
%----------------------------------
pt = get(get(obj,'Parent'),'CurrentPoint');
ud = get(gcf,'UserData');
ud.SIPos = round(pt(1,2));
ud.LRPos = round(pt(1,1));

set(ud.BackViewVC,'XData',[ud.LRPos ud.LRPos]);
set(ud.BackViewHC,'YData',[ud.SIPos ud.SIPos]);

set(ud.TopViewHandle,'CData',GetTopImage(ud,ud.SIPos));
set(ud.TopViewVC,'XData',[ud.LRPos ud.LRPos]);

set(ud.SideViewHandle,'CData',GetSideImage(ud,ud.LRPos));
set(ud.SideViewHC,'YData',[ud.SIPos ud.SIPos]);

UpdatePositionInfo(ud,handles);
ShowTalInfo(ud,handles);
set(gcf,'UserData',ud);

function TopViewKlicked(obj,event,handles)
%----------------------------------
pt = get(get(obj,'Parent'),'CurrentPoint');
ud = get(gcf,'UserData');
ud.APPos = round(pt(1,2));
ud.LRPos = round(pt(1,1));

set(ud.TopViewHC,'YData',[ud.APPos ud.APPos]);
set(ud.TopViewVC,'XData',[ud.LRPos ud.LRPos]);

set(ud.SideViewHandle,'CData',GetSideImage(ud,ud.LRPos));
set(ud.SideViewVC,'XData',[ud.APPos ud.APPos]);

set(ud.BackViewHandle,'CData',GetBackImage(ud,ud.APPos));
set(ud.BackViewVC,'XData',[ud.LRPos ud.LRPos]);

UpdatePositionInfo(ud,handles);
ShowTalInfo(ud,handles);
set(gcf,'UserData',ud);

function ShowColorScale(handles,ud)
subplot(handles.Colorbar);
y = 129:255;
cutoff = floor(ud.Threshold / ud.Max * 64);
y((64-cutoff):(64+cutoff)) = 0;
y = repmat(y',[1,10]);
imagesc(y,[0 255]);
%axis xy
axis(handles.Colorbar,'off','xy');

function editMax_Callback(hObject, eventdata, handles)
% hObject    handle to editMax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

ud = get(handles.figure1,'UserData');
ud.Max = str2double(get(handles.editMax,'String'));

set(ud.SideViewHandle,'CData',GetSideImage(ud,ud.LRPos));
set(ud.BackViewHandle,'CData',GetBackImage(ud,ud.APPos));
set(ud.TopViewHandle,'CData',GetTopImage(ud,ud.SIPos));
ShowColorScale(handles,ud);
set(handles.figure1,'UserData',ud);


% --- Executes during object creation, after setting all properties.
function editMax_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editMax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function editThreshold_Callback(hObject, eventdata, handles)
% hObject    handle to editThreshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

ud = get(handles.figure1,'UserData');
ud.Threshold = str2double(get(handles.editThreshold,'String'));

set(ud.SideViewHandle,'CData',GetSideImage(ud,ud.LRPos));
set(ud.BackViewHandle,'CData',GetBackImage(ud,ud.APPos));
set(ud.TopViewHandle,'CData',GetTopImage(ud,ud.SIPos));
ShowColorScale(handles,ud);
set(handles.figure1,'UserData',ud);


% --- Executes during object creation, after setting all properties.
function editThreshold_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editThreshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function EditXPos_Callback(hObject, eventdata, handles)
% hObject    handle to EditXPos (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

ud = get(handles.figure1,'UserData');
ud.LRPos = round(str2double(get(hObject,'String'))) + ud.XOff;
set(ud.SideViewHandle,'CData',GetSideImage(ud,ud.LRPos));
set(ud.BackViewVC,'XData',[ud.LRPos ud.LRPos]);
set(ud.TopViewVC,'XData',[ud.LRPos ud.LRPos]);
set(handles.figure1,'UserData',ud);


% --- Executes during object creation, after setting all properties.
function EditXPos_CreateFcn(hObject, eventdata, handles)
% hObject    handle to EditXPos (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function EditYPos_Callback(hObject, eventdata, handles)
% hObject    handle to EditYPos (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

ud = get(handles.figure1,'UserData');
ud.APPos = round(str2double(get(hObject,'String'))) + ud.YOff;
set(ud.BackViewHandle,'CData',GetBackImage(ud,ud.APPos));
set(ud.SideViewVC,'XData',[ud.APPos ud.APPos]);
set(ud.TopViewHC,'YData',[ud.APPos ud.APPos]);
set(handles.figure1,'UserData',ud);



% --- Executes during object creation, after setting all properties.
function EditYPos_CreateFcn(hObject, eventdata, handles)
% hObject    handle to EditYPos (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function EditZPos_Callback(hObject, eventdata, handles)
% hObject    handle to EditZPos (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

ud = get(handles.figure1,'UserData');
ud.SIPos = round(str2double(get(hObject,'String'))) + ud.ZOff;
set(ud.TopViewHandle,'CData',GetTopImage(ud,ud.SIPos));

set(ud.BackViewHC,'YData',[ud.SIPos ud.SIPos]);
set(ud.SideViewHC,'YData',[ud.SIPos ud.SIPos]);
set(handles.figure1,'UserData',ud);




% --- Executes during object creation, after setting all properties.
function EditZPos_CreateFcn(hObject, eventdata, handles)
% hObject    handle to EditZPos (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function UpdatePositionInfo(ud,handles)

set(handles.EditXPos,'String',sprintf('%i',ud.LRPos-ud.XOff));
set(handles.EditYPos,'String',sprintf('%i',ud.APPos-ud.YOff));
set(handles.EditZPos,'String',sprintf('%i',ud.SIPos-ud.ZOff));

function ShowTalInfo(ud,handles)

dltX = ud.VoxInfo{1} - (ud.LRPos-ud.XOff);
dltY = ud.VoxInfo{2} - (ud.APPos-ud.YOff); %+4);
dltZ = ud.VoxInfo{3} - (ud.SIPos-ud.ZOff); %-3);

dlt = dltX.*dltX + dltY.*dltY + dltZ.* dltZ;

[d,idx] = min(dlt);

Lobe = ud.VoxInfo{4};
Gyr  = ud.VoxInfo{5};
BA   = ud.VoxInfo{6};

set(handles.dspVoxelNumber,'String',sprintf('%i',idx));
Match = sprintf('Closest voxel is %2.1f mm away',sqrt(d));
StatInfo = sprintf('Value: %5.3f',ud.Lori(idx));
set(handles.TextInfo,'String',{Match StatInfo 'Location:' Lobe{idx},Gyr{idx},BA{idx}});


% --- Executes on button press in ButtonMax.
function ButtonMax_Callback(hObject, eventdata, handles)
% hObject    handle to ButtonMax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
DoMinMax(handles,1);

function DoMinMax(handles,flag)

ud = get(handles.figure1,'UserData');
if flag == 1
    [d,idx] = max(ud.Lori);
else
    [d,idx] = min(ud.Lori);
end

X = ud.VoxInfo{1};
Y = ud.VoxInfo{2};
Z = ud.VoxInfo{3};
ud.LRPos = X(idx)+ud.XOff;
ud.APPos = Y(idx)+ud.YOff; %-4;
ud.SIPos = Z(idx)+ud.ZOff; %+3;

set(ud.SideViewHandle,'CData',GetSideImage(ud,ud.LRPos));
set(ud.BackViewHandle,'CData',GetBackImage(ud,ud.APPos));
set(ud.TopViewHandle,'CData',GetTopImage(ud,ud.SIPos));

set(ud.BackViewHC,'YData',[ud.SIPos ud.SIPos]);
set(ud.BackViewVC,'XData',[ud.LRPos ud.LRPos]);

set(ud.SideViewHC,'YData',[ud.SIPos ud.SIPos]);
set(ud.SideViewVC,'XData',[ud.APPos ud.APPos]);

set(ud.TopViewHC,'YData',[ud.APPos ud.APPos]);
set(ud.TopViewVC,'XData',[ud.LRPos ud.LRPos]);

Lobe = ud.VoxInfo{4};
Gyr  = ud.VoxInfo{5};
BA   = ud.VoxInfo{6};

set(handles.dspVoxelNumber,'String',sprintf('%i',idx));
Match = sprintf('Closest voxel is %2.1f mm away',0);
StatInfo = sprintf('Value: %5.3f',ud.Lori(idx));
set(handles.TextInfo,'String',{Match StatInfo 'Location:' Lobe{idx},Gyr{idx},BA{idx}});
UpdatePositionInfo(ud,handles);

set(handles.figure1,'UserData',ud);

function TextInfo_Callback(hObject, eventdata, handles)
% hObject    handle to TextInfo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of TextInfo as text
%        str2double(get(hObject,'String')) returns contents of TextInfo as a double


% --- Executes during object creation, after setting all properties.
function TextInfo_CreateFcn(hObject, eventdata, handles)
% hObject    handle to TextInfo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in ButtonMin.
function ButtonMin_Callback(hObject, eventdata, handles)
% hObject    handle to ButtonMin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
DoMinMax(handles,0);

function h = LorBlueRed(m)

m = m - 2;
n = fix(1/2*m);

r = [zeros(m-n,1);[0;0]; (1:n)'/n    ];
g = [zeros(m-n,1);[0;0]; zeros(m-n,1)];
b = [(n:-1:1)'/n ;[0;0]; zeros(m-n,1)];    

h = [r g b];


% --- Executes on selection change in lbAnatomies.
function lbAnatomies_Callback(hObject, eventdata, handles)
% hObject    handle to lbAnatomies (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns lbAnatomies contents as cell array
%        contents{get(hObject,'Value')} returns selected item from lbAnatomies

ud = get(handles.figure1,'UserData');
idx = get(hObject,'Value');
path = get(hObject,'UserData');

fnames = get(hObject,'String');
fp = fopen([path '\' fnames{idx}]);
ud.Anatomy = reshape(fread(fp,'uint8'),[181,217,181]) / 2;
fclose(fp);

set(handles.figure1,'UserData',ud);

set(ud.SideViewHandle,'CData',GetSideImage(ud,ud.LRPos));
set(ud.BackViewHandle,'CData',GetBackImage(ud,ud.APPos));
set(ud.TopViewHandle,'CData',GetTopImage(ud,ud.SIPos));



% --- Executes during object creation, after setting all properties.
function lbAnatomies_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lbAnatomies (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in DumpToSLOR.
function DumpToSLOR_Callback(hObject, eventdata, handles)
% hObject    handle to DumpToSLOR (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[name,path] = uiputfile('*.slor','Save as SLOR file');
if name == 0
    return
end

ud = get(handles.figure1,'UserData');
fid = fopen(fullfile(path,name),'wb');
fullfile(path,name)
fwrite(fid,ud.Lori,'float32');
fclose(fid);



function dspVoxelNumber_Callback(hObject, eventdata, handles)
% hObject    handle to dspVoxelNumber (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of dspVoxelNumber as text
%        str2double(get(hObject,'String')) returns contents of dspVoxelNumber as a double


% --- Executes during object creation, after setting all properties.
function dspVoxelNumber_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dspVoxelNumber (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function SaveBitmap(v,map)

[filename, pathname] = uiputfile('image.bmp', 'Save image as');
imwrite(flipud(v),map,fullfile(pathname,filename),'bmp');


% --- Executes on button press in SaveBackView.
function SaveBackView_Callback(hObject, eventdata, handles)
% hObject    handle to SaveBackView (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

ud = get(handles.figure1,'UserData');
v = GetBackImage(ud,ud.APPos);
SaveBitmap(v,ud.colormap);

% --- Executes on button press in SaveSideView.
function SaveSideView_Callback(hObject, eventdata, handles)
% hObject    handle to SaveSideView (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

ud = get(handles.figure1,'UserData');
v = GetSideImage(ud,ud.LRPos);
FreezeColors();
get(ud.SideViewHandle)
SaveBitmap(v,ud.colormap);

% --- Executes on button press in SaveTopView.
function SaveTopView_Callback(hObject, eventdata, handles)
% hObject    handle to SaveTopView (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

ud = get(handles.figure1,'UserData');
v = GetTopImage(ud,ud.SIPos);
SaveBitmap(v,ud.colormap);
