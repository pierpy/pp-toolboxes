function varargout = Ragu(varargin)
% Ragu M-file for Ragu.fig
%      Ragu, by itself, creates a new Ragu or raises the existing
%      singleton*.
%
%      H = Ragu returns the handle to a new Ragu or the handle to
%      the existing singleton*.
%
%      Ragu('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in Ragu.M with the given input arguments.
%
%      Ragu('Property','Value',...) creates a new Ragu or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Ragu_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Ragu_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Ragu

% Copyright 2009-2011 Thomas Koenig
% distributed under the terms of the GNU AFFERO General Public License

% Last Modified by GUIDE v2.5 20-Jun-2016 13:09:54

% Begin initialization code - DO NOT EDIT

gui_Singleton = 0;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Ragu_OpeningFcn, ...
                   'gui_OutputFcn',  @Ragu_OutputFcn, ...
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

% --- Executes just before Ragu is made visible.
function Ragu_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Ragu (see VARARGIN)

% Choose default command line output for Ragu
handles.output = hObject;
handles.CurrentView = 'None';

% Update handles structure
guidata(hObject, handles);

if numel(varargin) > 0
    OpenMenuItem_Callback(hObject, eventdata, handles, varargin{1});
end

% This sets up the initial plot - only do when we are invisible
% so window can get raised using Ragu.

% UIWAIT makes Ragu wait for user response (see UIRESUME)
% uiwait(handles.figure1);
clc

if verLessThan('matlab', '7.10')
    rand('seed',sum(100*clock));
elseif verLessThan('matlab','7.14')
    s = RandStream.create('mt19937ar','seed',sum(100*clock));
    RandStream.setDefaultStream(s);
else
    s = RandStream.create('mt19937ar','seed',sum(100*clock));
    RandStream.setGlobalStream(s);

%data.CurrentView = 'None';
%guidata(hObject,data);    
    
end



% --- Outputs from this function are returned to the command line.
function varargout = Ragu_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --------------------------------------------------------------------
function FileMenu_Callback(hObject, eventdata, handles)
% hObject    handle to FileMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --------------------------------------------------------------------
function OpenMenuItem_Callback(hObject, eventdata, handles, fullfilename)
% hObject    handle to OpenMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
out = get(handles.output,'UserData');
CheckForResultsToSave(out);

if nargin < 4
    [filename, pathname] = uigetfile('*.mat', 'Load data from Matlab file');
    if isequal(filename,0) || isequal(pathname,0)
        return
    end
    cd(pathname);
    fullfilename = fullfile(pathname, filename);
end

if ~isempty(whos('-file',fullfilename,'RandomizerData'))
    load(fullfilename,'RandomizerData');
    rd  = RandomizerData;
elseif ~isempty(whos('-file',fullfilename,'rd'))
    load(fullfilename,'rd');
else
    uiwait(errordlg('No valid data in that file','Load data from Matlab file','modal'));
end

set(handles.figure1,'Name',['Ragu: ' fullfilename]);

if (~isfield(rd,'strF1'))
    rd.strF1 = 'Factor 1';
end

if (~isfield(rd,'strF2'))
    rd.strF2 = 'Factor 2';
end

rd = SetDefaults(rd);

if (~isfield(rd,'GFPPTanova'))
    rd.GFPPTanova = [];
end

if (~isfield(rd,'GFPTanovaEffectSize'))
    rd.GFPTanovaEffectSize = [];
end

if (~isfield(rd,'DLabels1'))
    uiwait(warndlg('Set level names before continuing','Load data from Matlab file','modal'));
    set(handles.output,'UserData',rd);
    MenuWithDesign_Callback(hObject, eventdata, handles);
    rd = get(handles.output,'UserData');
end

if (~isfield(rd,'GroupLabels'))
    i = unique(rd.IndFeature);
    if (numel(i) == 1 && rd.ContBetween == false)
        rd.GroupLabels = {'Group1'};
    else
        uiwait(warndlg('Set group names before continuing','Load data from Matlab file','modal'));
        set(handles.output,'UserData',rd);
        MenuBetweenDesign_Callback(hObject, eventdata, handles);
        rd = get(handles.output,'UserData');
    end
end



if (isfield(rd,'SamplingRate'))
    if rd.FreqDomain == 0
        rd.DeltaX = 1000 / rd.SamplingRate;
        rd.txtX = 'ms';
    else
        rd.DeltaX = rd.SamplingRate / 2 / size(rd.V,4);
        rd.txtX = 'Hz';
    end
    rd = rmfield(rd,'SamplingRate');
    rd = rmfield(rd,'FreqDomain');
end
    
if isfield(rd,'Modified')
    rd = rmfield(rd,'Modified');
end

if isfield(rd,'EndFrame')
    if rd.EndFrame > size(rd.V,4)
        rd.EndFrame = size(rd.V,4);
    end
end

if isfield(rd,'StartFrame')
    if rd.StartFrame > rd.EndFrame
        rd.StartFrame = rd.EndFrame;
    end
end

if size(rd.V,2) == 1
    rd.Design = [1 1];
end

rd.FileName = fullfilename;
set(handles.output,'UserData',rd);
cld();
Randomizer_ShowDataSet2();
UpdateGUI(handles);

% --------------------------------------------------------------------
function PrintMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to PrintMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

printdlg(handles.figure1);

% --------------------------------------------------------------------
function CloseMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to CloseMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
CheckForResultsToSave(get(handles.figure1,'UserData'));
delete(handles.figure1);

% --------------------------------------------------------------------
function SaveMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to SaveMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%UpdateParameters(handles);
rd = get(handles.output,'UserData');
[filename, pathname] = uiputfile('*.mat', 'Save data to Matlab file');
if isequal(filename,0) || isequal(pathname,0)
    return
end
if isfield(rd,'Modified')
    rd = rmfield(rd,'Modified');
end

fullfilename = fullfile(pathname, filename);
%rd.History = [rd.History,{['Data saved as: ' fullfilename]}];

save(fullfilename,'rd');

%cd(pathname);
rd.FileName = fullfilename;
set(handles.output,'UserData',rd);
set(handles.figure1,'Name',['Ragu: ' fullfilename]);

% --------------------------------------------------------------------
function MenuImport_Callback(hObject, eventdata, handles)
% hObject    handle to MenuImport (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
out = get(handles.output,'UserData');

if ~isempty(out)
    disp('hi');
    out = ResetResults(out);
    
    if isempty(out)
        return;
    end
else
    out = InitData();
end

Output = Randomizer_Import_ERP(hObject, eventdata, handles,out);
if isempty(Output)
     return;
end
x = get(Output,'UserData');
%out = get(handles.output,'UserData');
out.V = x.V;
out.Names = x.Names;
out.conds = x.conds;
out.Directory = x.Directory;
out.Mask      = x.Mask;
out.MakeAR    = x.MakeAR;

if (size(out.IndFeature,1) ~= size(out.Names,1))
    out.IndFeature = ones(size(out.Names,1),1);
    out.GroupLabels = {''};
end

if (isfield(out,'Channel'))
    if numel(out.Channel) ~= size(out.V,3)
        uiwait(warndlg('Loaded data is incompatible with the current montage, montage has been cleared'));
        out = rmfield(out,'Channel');
    end
end

if size(out.V,2) == 1
    out.Design = [1 1];
end

set(handles.output,'UserData',out);
close(Output);

%dummy = ['Data import from: ', out.Directory]; 
%out.History = [out.History,{dummy}];
%dummy = ['File mask: ', out.Mask]; 
%out.History = [out.History,{dummy}];
%dummy = 'Conditions: ';
%for c = 1:numel(out.conds)
%    dummy = [dummy out.conds{c} ' '];
%end

%out.History = [out.History,{dummy}];

cld();

UpdateGUI(handles);
MenuDataProperties_Callback([],[],handles);
MenuMontage_Callback([],[],handles);
Randomizer_ShowDataSet2();

% --------------------------------------------------------------------
function MenuMontage_Callback(hObject, eventdata, handles)
% hObject    handle to MenuMontage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%d = UpdateParameters(handles);
d = get(handles.output,'UserData');
[filename, pathname,FilterIndex] = uigetfile({'*.mat','Matlab structure from Analyzer (*.mat)';'*.txt','ASCII Text File with XYZ coordinates';'*.xyz','Cartool XYZ coordinates';'*.sxyz','sLORETA XYZ file'}, 'Load montage from Matlab file');
if isequal(filename,0) || isequal(pathname,0)
    return
end
switch FilterIndex
    case 1
        if isempty(whos('-file',fullfile(pathname, filename),'Channel'))
            uiwait(errordlg('No valid data in that file','Load montage from Matlab file','modal'));
        end
        load(fullfile(pathname, filename),'Channel');
    
        if numel(Channel) ~= size(d.V,3)
            uiwait(errordlg('Channel number mismatch','Load montage from Matlab file','modal'));
        else
            d.Channel = Channel;    
%            d.History = [d.History,{['Montage loaded from: ', fullfile(pathname, filename)]}];
            set(handles.output,'UserData',d);
        end
    case 2
        ChanPos = load(fullfile(pathname, filename));
        if (size(ChanPos,1) == 3 && size(ChanPos,2) == size(d.V,3))
            d.Channel = ChanPos;
            set(handles.output,'UserData',d);
        elseif (size(ChanPos,2) == 3 && size(ChanPos,1) == size(d.V,3))
                d.Channel = ChanPos';
%                d.History = [d.History,{['Montage loaded from: ', fullfile(pathname, filename)]}];

                set(handles.output,'UserData',d);
        else
            uiwait(errordlg('Channel number mismatch','Load montage from Matlab file','modal'));
        end
    case 3    
        ChanPos = ReadXYZ(fullfile(pathname, filename));
        if ~isempty(ChanPos)
            d.Channel = ChanPos';
%            d.History = [d.History,{['Montage loaded from XYZ: ', fullfile(pathname, filename)]}];

        end
        set(handles.output,'UserData',d);
        
    case 4    
        ChanPos = ReadSXYZ(fullfile(pathname, filename));
        if ~isempty(ChanPos)
            d.Channel = ChanPos;
%            d.History = [d.History,{['Montage loaded from SXYZ: ', fullfile(pathname, filename)]}];

        end
        set(handles.output,'UserData',d);

end

cld();

% --------------------------------------------------------------------
function MenuDesign_Callback(hObject, eventdata, handles)
% hObject    handle to MenuDesign (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --------------------------------------------------------------------
function MenuWithDesign_Callback(hObject, eventdata, handles)
% hObject    handle to MenuWithDesign (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
out = get(handles.output,'UserData');

if(isempty(out.V))
    uiwait(errordlg('Load subject data first','Define within subject design','modal'));
    return;
end

out = ResetResults(out);

if isempty(out)
    return
end

Output = Randomizer_Design(hObject, eventdata, handles,out);
if isempty(Output)
    return;
end

out = get(Output,'UserData');
cld();
%out.History = [out.History,{'Within subject design changed.'}];

set(handles.output,'UserData',out);
close(Output);
UpdateGUI(handles);

% --------------------------------------------------------------------
function MenuBetweenDesign_Callback(hObject, eventdata, handles)
% hObject    handle to MenuBetweenDesign (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
out = get(handles.output,'UserData');

if(isempty(out.V))
    uiwait(errordlg('Load ERP/EEG data first','Define behavioral data / groups','modal'));
    return;
end

out = ResetResults(out);

if isempty(out)
    return;
end

Output = Randomizer_IndFeatures(hObject, eventdata, handles,out);

if isempty(Output)
    return;
end
out = get(Output,'UserData');
%out.History = [out.History,{'Within subject design changed.'}];

set(handles.output,'UserData',out);
close(Output);
UpdateGUI(handles);

% --------------------------------------------------------------------
function MenuData_Callback(hObject, eventdata, handles)
% hObject    handle to MenuData (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function MenuAnalysis_Callback(hObject, eventdata, handles)
% hObject    handle to MenuAnalysis (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function MenuTCT_Callback(hObject, eventdata, handles)
% hObject    handle to MenuTCT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
out = get(handles.figure1,'UserData');


Output = Randomizer_Tanova_Options(hObject, eventdata, handles,out,false);
out = get(Output,'UserData');

if ~isfield(out,'Continue')
    return
end

close(Output);

out = rmfield(out,'Continue');

if TestForSingleGroups(handles)
    return;
end

out = ComputeAnTopography(out);
out.Modified = true;
%out = CheckForResultsToSave(out);

%out.History = [out.History,{sprintf('TCT computed from %i to %i frames',out.StartFrame,out.EndFrame)}];

set(handles.output,'UserData',out);

cld();
ShowAnTopographyResults(out,handles.figure1);

% --------------------------------------------------------------------



function MenuTanova_Callback(hObject, eventdata, handles)
% hObject    handle to MenuTanova (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

DoTheTanova(hObject, eventdata, handles,false);

function DoTheTanova(hObject, eventdata, handles,DoTheGFP)
out = get(handles.output,'UserData');
% UpdateParameters(handles);


if TestForSingleGroups(handles)
    return;
end

Output = Randomizer_Tanova_Options(hObject, eventdata, handles,out,false);
out = get(Output,'UserData');

if ~isfield(out,'Continue')
    return
end

out = rmfield(out,'Continue');

if isfield(out,'CritFDR_p')
    out = rmfield(out,'CritFDR_p');
end

set(handles.output,'UserData',out);
close(Output);

if DoTheGFP == false
        out.DoGFP = 0;
else
        out.DoGFP = 1;
end

out = Randomizer_ComputeTanova(out);

if ~isstruct(out)
    return
end

if DoTheGFP == false
    out.stt = RaguSingleThresholdTest(out);
    out.CritDuration = ones(2,4) * 1000000000000000;
    out.TanovaHits = [];
    out.PHitCount = [];
else
    out.CritDurationGFP = ones(2,4) * 1000000000000000;
    out.GFPHits = [];
    out.PHitCountGFP = [];
end
out.Modified = true;
set(handles.output,'UserData',out);
%CheckForResultsToSave(out);
cld();
UpdateGUI(handles);
Randomizer_ShowTanovaResults(out,handles.figure1);

% --------------------------------------------------------------------
function MenuCountStats_Callback(hObject, eventdata, handles)
% hObject    handle to MenuCountStats (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

d = get(handles.output,'UserData');

d.TanovaHits = squeeze(sum(d.PTanova(:,:,d.StartFrame:d.EndFrame,:) < d.Threshold,3));
for i = 1:2
    for j = 1:4
        d.PHitCount(i,j) = squeeze(sum(d.TanovaHits(i,j,:)>= d.TanovaHits(i,j,1),3)) / size(d.TanovaHits,3);
    end
end

set(handles.output,'UserData',d);
Randomizer_ShowHitCountResults(d);

% --------------------------------------------------------------------
function MenuCountStatsGFP_Callback(hObject, eventdata, handles)
% hObject    handle to MenuCountStatsGFP (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

d = get(handles.output,'UserData');

d.GFPHits = squeeze(sum(d.GFPPTanova(:,:,d.StartFrame:d.EndFrame,:) < d.Threshold,3));
for i = 1:2
    for j = 1:4
        d.PHitCountGFP(i,j) = squeeze(sum(d.GFPHits(i,j,:)>= d.GFPHits(i,j,1),3)) / size(d.GFPHits,3);
    end
end

set(handles.output,'UserData',d);
Randomizer_ShowHitCountResults(d,3);


% --------------------------------------------------------------------
function MenuClustStats_Callback(hObject, eventdata, handles)
% hObject    handle to MenuClustStats (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

d = get(handles.output,'UserData');
hit = d.PTanova(:,:,d.StartFrame:d.EndFrame,:) < d.Threshold;
for i = 1:2
    for j = 1:4
        d.PHitDuration{i,j} = RaguClustSizeV2(squeeze(hit(i,j,:,:)));
        p = 1-(cumsum(d.PHitDuration{i,j})/sum(d.PHitDuration{i,j}));
        crit = find(p <= d.Threshold);
        if ~isempty(crit)
            d.CritDuration(i,j) = crit(1);
        else
            d.CritDuration(i,j) = size(d.V,4);
        end
    end
end

set(handles.output,'UserData',d);
TanovaResults_Callback(hObject, eventdata, handles);

Randomizer_ShowHitDurationResults(d);

% --------------------------------------------------------------------
function MenuClustStatsGFP_Callback(hObject, eventdata, handles)
% hObject    handle to MenuClustStatsGFP (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

d = get(handles.output,'UserData');
hit = d.GFPPTanova(:,:,d.StartFrame:d.EndFrame,:) < d.Threshold;
for i = 1:2
    for j = 1:4
        d.PHitDurationGFP{i,j} = RaguClustSizeV2(squeeze(hit(i,j,:,:)));
        p = 1-(cumsum(d.PHitDurationGFP{i,j})/sum(d.PHitDurationGFP{i,j}));
        crit = find(p <= d.Threshold);
        if ~isempty(crit)
            d.CritDurationGFP(i,j) = crit(1);
        else
            d.CritDurationGFP(i,j) = size(d.V,4);
        end
    end
end

set(handles.output,'UserData',d);
View_GFP_Callback(hObject, eventdata, handles);

Randomizer_ShowHitDurationResults(d,1);






% --------------------------------------------------------------------
function MenuRNDOptions_Callback(hObject, eventdata, handles)
% hObject    handle to MenuRNDOptions (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

dat = get(handles.output,'UserData');

h = RaguRNDOpt2(hObject, eventdata, handles,dat);

if isempty(h)
    return
end

out = get(h,'UserData');

dat.Threshold = out.Threshold;
dat.NoXing = out.NoXing;
dat.Iterations = out.Iterations;
dat.Normalize = out.Normalize;
close(h);

set(handles.output,'UserData',dat);

% --------------------------------------------------------------------
function MenuDataProperties_Callback(hObject, eventdata, handles)
% hObject    handle to MenuDataProperties (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
dat = get(handles.output,'UserData');

h = RaguDataProps(hObject, eventdata, handles,dat);

if (isempty(h))
    return
end

out = get(h,'UserData');
dat.TimeOnset = out.TimeOnset;
dat.DeltaX = out.DeltaX;
dat.txtX = out.txtX;

close(h);
cld();
set(handles.output,'UserData',dat);

% --------------------------------------------------------------------
function MenuViewData_Callback(hObject, eventdata, handles)
% hObject    handle to MenuViewData (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%out = UpdateParameters(handles);
%Output = Randomizer_ShowDataSet2(hObject, eventdata, handles,out);
cld();
Randomizer_ShowDataSet2();

% --------------------------------------------------------------------
function RandomizerData = CheckForResultsToSave(out)

RandomizerData = out;
if ~isfield(out,'Modified')
    return;
end

resp  = questdlg('What would you like to do?','Unsaved data in memory','Save results','Clear results','Cancel','Save results');
    
switch resp
    case 'Save results'
        
        [filename, pathname] = uiputfile('*.mat', 'Save data to Matlab file');
        if filename ~= 0
            rd = rmfield(RandomizerData,'Modified');
            save(fullfile(pathname, filename),'rd');
            cd(pathname);
        end
    case 'Clear results'
        return
        
    case 'Cancel'
        rd = [];
        return
end

% --------------------------------------------------------------------

function out = InitData(in)

if nargin > 0
    out = in;
end

out.V = [];
out.IndFeature = [];
out.Design = [];
%out.History = {};
out = InitResults(out);


function out = InitResults(in)

if nargin > 0
    out = in;
end
out.PTanova = [];
out.TanovaEffectSize = [];
out.GFPPTanova = [];
out.GFPTanovaEffectSize = [];
out.MeanGFP = [];
out.pTopCons = [];
out.TanovaHits = [];
out.CritDuration = ones(2,4) * 1000000000000000;
out.PHitCount = [];
out.PHitDuration = [];
out.MSMaps = [];
out.MSEffectSize = [];
out.pMSStats = [];

out = SetDefaults(out);


function out = SetDefaults(out)
out.Threshold          = DefaultTo(out,'Threshold',.05);
out.Iterations         = DefaultTo(out,'Iterations',5000);
out.strF1              = DefaultTo(out,'strF1','Factor 1');
out.strF2              = DefaultTo(out,'strF2','Factor 2');
out.IndName            = DefaultTo(out,'IndName','Group');
out.MeanInterval       = DefaultTo(out,'MeanInterval',0);
out.Normalize          = DefaultTo(out,'Normalize',1);
out.FrequencyDomain    = DefaultTo(out,'FrequencyDomain',0);
out.DoAnTopFFT         = DefaultTo(out,'DoAnTopFFT',0);
out.NoXing             = DefaultTo(out,'NoXing',false);
out.txtX               = DefaultTo(out,'txtX','ms');
out.DeltaX             = DefaultTo(out,'DeltaX',1);
out.TimeOnset          = DefaultTo(out,'TimeOnset',0);
out.StartFrame         = DefaultTo(out,'StartFrame',[]);
out.EndFrame           = DefaultTo(out,'EndFrame',[]);
out.ContBetween        = DefaultTo(out,'ContBetween',false);
out.BarGraph           = DefaultTo(out,'BarGraph',0);
out.TwoFactors         = DefaultTo(out,'TwoFactors',0);
out.ContF1             = DefaultTo(out,'ContF1',0);
out.MapStyle           = DefaultTo(out,'MapStyle',2);
out.DoGFP              = DefaultTo(out,'DoGFP',0);
out.OptStart           = DefaultTo(out,'OptStart',3);
out.OptEnd             = DefaultTo(out,'OptEnd',10);
out.OptNTraining       = DefaultTo(out,'OptNTraining',floor(sum(~isnan(out.IndFeature)))/2);
out.FixedK             = DefaultTo(out,'FixedK',0);
out.NoInconsistentMaps = DefaultTo(out,'NoInconsistentMaps',0);
out.FDR                = DefaultTo(out,'FDR',0.2);
out.XValRestarts       = DefaultTo(out,'XValRestarts',50);


function res = DefaultTo(in,field,val)

if(isfield(in,field))
    res = getfield(in,field);
else
    res = val;
end


% --------------------------------------------------------------------
function r = ResetResults(out)

out = CheckForResultsToSave(out);

if isempty(out)
    r = [];
    return
end
cld();

if isfield(out,'Modified')
    rmfield(out,'Modified');
end

out.PTanova = [];
out.pMSStats = [];
out.TanovaEffectSize = [];
out.GFPPTanova = [];
out.GFPTanovaEffectSize = [];
out.MeanGFP = [];
out.pTopCons = [];
out.PHitDuration = [];
out.TanovaHits = [];
out.PHitCount = [];
out.CritDuration = ones(2,4) * 1000000000000000;

r = out;

% --------------------------------------------------------------------
function View_Callback(hObject, eventdata, handles)
% hObject    handle to View (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --------------------------------------------------------------------
function TanovaResults_Callback(hObject, eventdata, handles)
% hObject    handle to TanovaResults (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

out = get(handles.output,'UserData');
out.DoGFP = 0;
set(handles.output,'UserData',out);
if ~isempty(out.PTanova)
    cld();

    out.TanovaMeanEffectSize = squeeze(mean(out.TanovaEffectSize,3));
    Rank = 1:size(out.TanovaEffectSize,4);
    for i = 1:2
        for j = 1:4
            [mx,order] = sort(squeeze(out.TanovaMeanEffectSize(i,j,:)),'descend');
            out.PTanovaOverall(i,j,order) = Rank;
        end
    end
    
    out.PTanovaOverall = out.PTanovaOverall / out.Iterations;
    
    Randomizer_ShowTanovaResults(out,handles.figure1);
else
    switch questdlg('Run the analysis now?','Results not yet computed','Yes','No','Yes');
        case 'Yes'
            MenuTanova_Callback(hObject, eventdata, handles)
    end
end

% --------------------------------------------------------------------
function cld()

dt = findobj('Tag','DataTable');

if ~isempty(dt)
    delete(dt);
end


dt = findall(gcf,'Tag','AnTopTextBox');

if ~isempty(dt)
    delete(dt);
end

dt = findall(gcf,'Tag','Transient');

if ~isempty(dt)
    delete(dt);
end

dt = findall(gcf,'Tag','MSStatsTable');

if ~isempty(dt)
    delete(dt);
end



h = subplot('Position',[0 0 1 1]);
cla(h);
delete(h);

% --------------------------------------------------------------------
function UpdateGUI(handles)

out = get(handles.figure1,'UserData');
if(~isempty(out.V))
    DataThere = 'on';
else
    DataThere = 'off';
end

if(~isempty(out.V) && ~isempty(out.Design))
    DDThere = 'on';
else
    DDThere = 'off';
end

if ~isempty(out.Design)
    if (size(out.V,2) > 1 && (var(out.Design(:,1)) == 0 ))
        DDThere = 'off';
    end
end

MSOK = 'off';

if strcmp(DDThere,'on')   & out.ContBetween == false & out.ContF1 == false
    MSOK = 'on';
end
    
if isempty(out.PTanova)
    ARthere = 'off';
else
    ARthere = 'on';
end

if isempty(out.GFPPTanova)
    GFPthere = 'off';
else
    GFPthere = 'on';
end




MSMaps = 'off';
if isfield(out,'MSMaps')
    if ~isempty(out.MSMaps)
        MSMaps = 'on';
    end
end

if out.MapStyle == 1
	set(handles.MenuViewPrettyMaps,'Checked','off');
else
    set(handles.MenuViewPrettyMaps,'Checked','on');
end

set(handles.Tools_RemoveBL,'Enable',DataThere);
set(handles.ToolsFilter,'Enable',DataThere);
set(handles.T_Mapper,'Enable',DataThere);
set(handles.MenuTanova,'Enable',DDThere);
set(handles.MenuMicrostates,'Enable',DDThere);
set(handles.MenuAnalysis_GFP,'Enable',DDThere);
set(handles.TanovaResults,'Enable',DDThere);
set(handles.View_GFP,'Enable',DDThere);
set(handles.MenuTCT,'Enable',DataThere);
set(handles.FileSave,'Enable',DataThere);
set(handles.PrintMenuItem,'Enable',DataThere);
set(handles.Tools_FindOutlier,'Enable',DataThere);

set(handles.SaveMeta,'Enable',DataThere);
set(handles.SaveBitmap,'Enable',DataThere);
set(handles.PrintMenuItem,'Enable',DataThere);
set(handles.FileSaveFigure,'Enable',DataThere);
set(handles.EditCopy,'Enable',DataThere);
set(handles.EditCopyBitmap,'Enable',DataThere);
set(handles.Edit_AddLabel,'Enable',DataThere);

set(handles.MenuDataProperties,'Enable',DataThere);
set(handles.MenuRNDOptions,'Enable',DataThere);
set(handles.MenuMontage,'Enable',DataThere);
set(handles.TCTResults,'Enable',DataThere);
set(handles.Data_Export,'Enable',DataThere);
set(handles.MenuWithDesign,'Enable',DataThere);
set(handles.MenuBetweenDesign,'Enable',DataThere);
set(handles.SaveMenuItem,'Enable',DataThere);
set(handles.MenuViewData,'Enable',DataThere);
set(handles.MS_Cormat,'Enable',DataThere);


set(handles.Tanova_Overall,'Enable',ARthere);
set(handles.MenuCountStats,'Enable',ARthere);
set(handles.MenuAUCStats,'Enable',ARthere);
set(handles.MenuClustStats,'Enable',ARthere);
set(handles.MenuFileExpTanova,'Enable',ARthere);
set(handles.MenuFileExpGFP,'Enable',GFPthere);


set(handles.MSCompute,'Enable',MSOK);
set(handles.MSLoadMaps,'Enable',MSOK);
set(handles.MS_Stats,'Enable',MSMaps);
set(handles.MSViewStats,'Enable',MSMaps);
set(handles.MS_SaveMaps,'Enable',MSMaps);

set(handles.GFP_Overall,'Enable',GFPthere);
set(handles.MenuCountStatsGFP,'Enable',GFPthere);
set(handles.MenuClustStatsGFP,'Enable',GFPthere);
set(handles.MenuAUCStatsGFP,'Enable',GFPthere);


if(isfield(out,'OrgData'))
    set(handles.ToolsUndo,'Enable','On');
else
    set(handles.ToolsUndo,'Enable','Off');
end

% --------------------------------------------------------------------
function TCTResults_Callback(hObject, eventdata, handles)
% hObject    handle to TCTResults (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
out = get(handles.figure1,'UserData');

if ~isempty(out.pTopCons)
    cld();
    ShowAnTopographyResults(out,handles.figure1);
else
    switch questdlg('Run the analysis now?','Results not yet computed','Yes','No','Yes');
        case 'Yes'
            MenuTCT_Callback(hObject, eventdata, handles);
    end
end   


% --------------------------------------------------------------------
function MenuEdit_Callback(hObject, eventdata, handles)
% hObject    handle to MenuEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function EditCopy_Callback(hObject, eventdata, handles)
% hObject    handle to EditCopy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
SaveTheFigure('-dmeta',[],[]);

% --------------------------------------------------------------------
function SaveMeta_Callback(hObject, eventdata, handles)
% hObject    handle to SaveMeta (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
SaveTheFigure('-dmeta','*.wmf', 'Save figure to Metafile');


function SaveTheFigure(Device,mask,comment)

if (isempty(mask))
    print(Device);
else
    [filename, pathname] = uiputfile(mask,comment);
    if isequal(filename,0) || isequal(pathname,0)
        return
    end
    print(Device,fullfile(pathname, filename));
end


% --------------------------------------------------------------------
function SaveBitmap_Callback(hObject, eventdata, handles)
% hObject    handle to SaveBitmap (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
SaveTheFigure('-dbitmap','*.bmp', 'Save figure to a bitmap');


% --------------------------------------------------------------------
function EditCopyBitmap_Callback(hObject, eventdata, handles)
% hObject    handle to EditCopyBitmap (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
SaveTheFigure('-dbitmap',[],[]);


% --------------------------------------------------------------------
function MenuFileExpTanova_Callback(hObject, eventdata, handles, GFPFlag)
% hObject    handle to MenuFileExpTanova (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


DoGFP = 0;
if nargin > 3
    DoGFP = GFPFlag;
end

[filename, pathname] = uiputfile('*.txt', 'Save as textfile');
if isequal(filename,0) || isequal(pathname,0)
    return
end

[fp,err] = fopen(fullfile(pathname, filename),'wt');

if fp == -1
    errordlg(err);
    return
end
    
d = get(handles.figure1,'UserData');

fprintf(fp,'Time');

for i = 1:d.ng
    for j = 1:d.nc
        if (i * j > 1)
            fprintf(fp,'\t%s',d.titles{d.ng*j+i-2});
        end
    end
end

if DoGFP == 0
    out = d.PTanova;
else
    out = d.GFPPTanova;
end

for t = 1:size(out,3);
    fprintf(fp,'\n%5.2f',d.time(t));
    for i = 1:d.ng
        for j = 1:d.nc
            if (i * j > 1)
                fprintf(fp,'\t%4.4f',out(i,j,t,1));
            end
        end
    end
end
    
fclose(fp);


% --------------------------------------------------------------------
function MenuHelp_Callback(hObject, eventdata, handles)
% hObject    handle to MenuHelp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



% --------------------------------------------------------------------
function HelpContent_Callback(hObject, eventdata, handles)
% hObject    handle to HelpContent (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function HelpAbout_Callback(hObject, eventdata, handles)
% hObject    handle to HelpAbout (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

helpdlg(RaguDate,'Version info');


% --- Executes during object creation, after setting all properties.
function figure1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object deletion, before destroying properties.
function figure1_DeleteFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
out = get(handles.figure1,'UserData');
CheckForResultsToSave(out);


% --------------------------------------------------------------------
function ToolbarLoad_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to ToolbarLoad (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
OpenMenuItem_Callback(hObject, eventdata, handles);


% --------------------------------------------------------------------
function ToolbarSave_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to ToolbarSave (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
SaveMenuItem_Callback(hObject, eventdata, handles);


% --------------------------------------------------------------------
function Toolbar_Print_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to Toolbar_Print (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
printdlg();


% --------------------------------------------------------------------
function FileSaveFigure_Callback(hObject, eventdata, handles)
% hObject    handle to FileSaveFigure (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[filename, pathname] = uiputfile('*.fig', 'Save figure');
if isequal(filename,0) || isequal(pathname,0)
    return
end

saveas(handles.figure1,fullfile(pathname, filename));
%hgsave(handles.figure1,fullfile(pathname, filename),'-v6');


% --------------------------------------------------------------------
function FileOpenFigure_Callback(hObject, eventdata, handles)
% hObject    handle to FileOpenFigure (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[filename, pathname] = uigetfile('*.fig', 'Save figure');

if isequal(filename,0) || isequal(pathname,0)
    return
end

h = openfig(fullfile(pathname, filename));
set(h,'Name',filename);


% --------------------------------------------------------------------
function DataClear_Callback(hObject, eventdata, handles)
% hObject    handle to DataClear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
out = get(handles.output,'UserData');
out = InitData(out);
cld();
set(handles.output,'UserData',out);


% --------------------------------------------------------------------
function STT_Callback(hObject, eventdata, handles)
% hObject    handle to STT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

out = get(handles.figure1,'UserData');
[junk,out.stt] = RaguSingleThresholdTest(out);
RaguSTT(hObject, eventdata, handles,out);


% --------------------------------------------------------------------
function MenuMicrostates_Callback(hObject, eventdata, handles)
% hObject    handle to MenuMicrostates (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --------------------------------------------------------------------
function MSLoadMaps_Callback(hObject, eventdata, handles)
% hObject    handle to MSLoadMaps (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
out = get(handles.figure1,'UserData');

[filename, pathname] = uigetfile('*.asc;*.txt', 'Load data from ASCII file');
if isequal(filename,0) || isequal(pathname,0)
    return
end

MSMaps = load(fullfile(pathname, filename));

if ~isempty(MSMaps)
    if size(MSMaps,2) ~= size(out.V)
        uiwait(errordlg('Channel number mismatch','Load data from Matlab file','modal'));
    else
        out.MSMaps = MSMaps;
        out.nStates = size(MSMaps,1);
   
        h = RaguMSParameters(hObject, eventdata, handles,out,1);
        if (isempty(h))
            return
        end
        
        u = get(h,'UserData');
        close(h);

        out.nReiter = u.nReiter;
        out.bSmoothLabels = u.bSmoothLabels;
        out.nWindowSize = u.nWindowSize;
        out.LabelPenalty = u.LabelPenalty;
        out.NoInconsistentMaps = u.NoInconsistentMaps;
  
        set(handles.figure1,'UserData',out);
        UpdateGUI(handles);
        cld();
        RaguShowMSFit(handles.figure1);
    end
end

% --------------------------------------------------------------------
function MS_Fit_Callback(hObject, eventdata, handles)
% hObject    handle to MS_Fit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



% --------------------------------------------------------------------
function MS_Stats_Callback(hObject, eventdata, handles)
% hObject    handle to MS_Stats (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
out = get(handles.figure1,'UserData');


Output = Randomizer_Tanova_Options(hObject, eventdata, handles,out,true);
out = get(Output,'UserData');

if ~isfield(out,'Continue')
    return
end

out = rmfield(out,'Continue');

set(handles.output,'UserData',out);
close(Output);

out = RaguMSStats(out);

set(handles.figure1,'UserData',out);
cld();
RaguShowMSFit(handles.figure1);



function maps = QuickComputeMSMaps(out,Subset,ProgBar,n_mods)

IndFeatureToUse = out.IndFeature;
IndFeatureToUse(Subset == 0) = NaN;
gmall2 = RaguGrandMeans(out.V(:,:,:,out.StartFrame:out.EndFrame),IndFeatureToUse);
   
gm = permute(gmall2,[1,2,4,3]);
gm = reshape(gm,[size(gm,1)*size(gm,2)*size(gm,3),size(gm,4)]);

if (out.NoInconsistentMaps == 1)
    if isempty(out.pTopCons)
        uiwait(errordlg('No valid TCT data available','Compute microstate maps','modal'));
        maps = NaN;
        return;
    end
    
    p = out.pTopCons(:,:,out.StartFrame:out.EndFrame);
    p = reshape(p,[size(p,1)*size(p,2)*size(p,3),1]);
    idx = p <= out.Threshold;
    gm = gm(idx,:);
end

if(out.UseAAHC)
    maps = RaguEEG_Mod_AAHC(gm,n_mods,ProgBar);
else
    for i = 1:numel(n_mods)
        maps{i} = NormDimL2(RaguEEG_Mod_r(gm,n_mods(i),out.nReiter,[],'p'),2);
    end
end


% --------------------------------------------------------------------
function out = ComputeMSMaps(out,Subset,ProgBar)

if isempty(Subset)
    gmall  = RaguGrandMeans(out.V(:,:,:,:),out.IndFeature);
    gmall2 = RaguGrandMeans(out.V(:,:,:,out.StartFrame:out.EndFrame),out.IndFeature);
else
    IndFeatureToUse = out.IndFeature;
    IndFeatureToUse(Subset == 0) = NaN;
    gmall2 = RaguGrandMeans(out.V(:,:,:,out.StartFrame:out.EndFrame),IndFeatureToUse);
    gmall = RaguGrandMeans(out.V(:,:,:,:),IndFeatureToUse);
end
    
gm = permute(gmall2,[1,2,4,3]);
gm = reshape(gm,[size(gm,1)*size(gm,2)*size(gm,3),size(gm,4)]);

if (out.NoInconsistentMaps == 1 && ~isempty(out.pTopCons))
    p = out.pTopCons(:,:,out.StartFrame:out.EndFrame);
    p = reshape(p,[size(p,1)*size(p,2)*size(p,3),1]);
    idx = p <= out.Threshold;
    gm = gm(idx,:);
end

if(out.UseAAHC)
    b_model = RaguEEG_Mod_AAHC(gm,out.nStates,ProgBar);
else
    [b_model] = RaguEEG_Mod_r(gm,out.nStates,out.nReiter,[],'p');
end

[MSClass,MSFit] = RaguFitMicrostates(gmall,b_model,out.bSmoothLabels,out.nWindowSize,out.LabelPenalty,out.StartFrame,out.EndFrame);

idx1 = 1:size(MSClass,1);
idx2 = 1:size(MSClass,2);
idx3 = 1:size(b_model,1);
on(idx1,idx2,idx3) = RaguMSOnOffsets(MSClass,size(b_model,1),MSFit,[], 0);

onr = reshape(on,size(on,1)*size(on,2),size(on,3));
m_on = zeros(size(onr,2),1);

for i = 1:size(onr,2)
    j = onr(:,i);
    j(isnan(j)) = [];
    m_on(i) = mean(j);
end

[s,ord] = sort(m_on);

out.MSMaps = b_model(ord,:);

if (isfield(out,'pMSStats'))
    out = rmfield(out,'pMSStats');
end



% --------------------------------------------------------------------
function MSCompute_Callback(hObject, eventdata, handles)
% hObject    handle to MSCompute (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

out = get(handles.figure1,'UserData');

Output = Randomizer_Tanova_Options(hObject, eventdata, handles,out,true);
out = get(Output,'UserData');
if isempty(out)
    return;
end
close(Output);

h = RaguMSParameters(hObject, eventdata, handles,out);
if (isempty(h))
    return
end

u = get(h,'UserData');
close(h);

out.nStates = u.nStates;
out.nReiter = u.nReiter;
out.bSmoothLabels = u.bSmoothLabels;
out.nWindowSize = u.nWindowSize;
out.LabelPenalty = u.LabelPenalty;
out.NoInconsistentMaps = u.NoInconsistentMaps;
out.FixedK   = u.FixedK;
out.OptStart = u.OptStart;
out.OptEnd = u.OptEnd;
out.XValRestarts = u.XValRestarts;
out.OptNTraining = u.OptNTraining;
out.UseAAHC = u.UseAAHC;
set(handles.figure1,'UserData',out);

if (out.FixedK == true)
    out = ComputeMSMaps(out,[],true);
else
    out = XValMSNumber(out);
end

if ~isstruct(out)
    return;
end

set(handles.figure1,'UserData',out);
figure(handles.figure1);

UpdateGUI(handles);
cld();
RaguShowMSFit(handles.figure1);


% --------------------------------------------------------------------
function MSViewStats_Callback(hObject, eventdata, handles)
% hObject    handle to MSViewStats (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

cld();
RaguShowMSFit(handles.figure1);


% --------------------------------------------------------------------
function MenuViewPrettyMaps_Callback(hObject, eventdata, handles)
% hObject    handle to MenuViewPrettyMaps (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

out = get(handles.figure1,'UserData');
if strcmp(get(hObject,'Checked'),'on')
    set(hObject,'Checked','off');
    out.MapStyle = 1;
else
    set(hObject,'Checked','on');
    out.MapStyle = 2;
end
set(handles.figure1,'UserData',out);


% --------------------------------------------------------------------
function MenuOptions_Callback(hObject, eventdata, handles)
% hObject    handle to MenuOptions (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Options_CLC_Callback(hObject, eventdata, handles)
% hObject    handle to Options_CLC (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
clc;


% --------------------------------------------------------------------
function ret = XValMSNumber(out)

nGroups = numel(unique(out.IndFeature(~isnan(out.IndFeature))));

res = nan(out.XValRestarts,out.OptEnd);
lres = nan(out.XValRestarts,out.OptEnd);
%MeanFit = nan(out.OptEnd);

if (out.NoInconsistentMaps == 1 && ~isempty(out.pTopCons))
    p = out.pTopCons;
    UsemeIdx = (p <= out.Threshold);
else
    UsemeIdx = ones(nGroups,size(out.V,2),size(out.V,4));
end

for run = 1:out.XValRestarts
    ReTry = 1;
    while ReTry == 1
        idx = randperm(numel(out.IndFeature));
        idx(isnan(out.IndFeature)) = nan;
        rnk = tiedrank(idx);
        
%        Grouping = randn(size(out.IndFeature));
        LearnGroup   = rnk <=out.OptNTraining;
        nGroupsLearn = numel(unique(out.IndFeature(LearnGroup)));
        TestGroup    = rnk  > out.OptNTraining;
        nGroupsTest  = numel(unique(out.IndFeature(TestGroup)));
        if (nGroupsLearn == nGroups) && (nGroupsTest == nGroups)
            ReTry = 0;
        end
    end
    IndFeatureToUse = out.IndFeature;
    IndFeatureToUse(LearnGroup == 0) = NaN;
    gmallL = RaguGrandMeans(out.V,IndFeatureToUse);

    IndFeatureToUse = out.IndFeature;
    IndFeatureToUse(TestGroup == 0) = NaN;
    gmallT = RaguGrandMeans(out.V,IndFeatureToUse);
  
    maps = QuickComputeMSMaps(out,LearnGroup,true,out.OptStart:out.OptEnd);
    for n = out.OptStart:out.OptEnd
        out.nStates = n;
        MapToUse = maps{n-out.OptStart+1};
        MSClass = RaguFitMicrostates(gmallL,MapToUse,out.bSmoothLabels,out.nWindowSize,out.LabelPenalty,out.StartFrame,out.EndFrame);

        [s1,s2,s3] = size(MSClass);
        gmallTPR = reshape(permute(gmallT,[1 2 4 3]),[s1*s2*s3,size(gmallT,3)]);
        gmallLPR = reshape(permute(gmallL,[1 2 4 3]),[s1*s2*s3,size(gmallL,3)]);
        MSClassR = reshape(MSClass,[s1*s2*s3,1]);
        model = zeros(numel(MSClassR),size(gmallL,3));
        UseMeR = reshape(UsemeIdx,[s1*s2*s3,1]);
        IsGood = ((MSClassR ~= 0) & UseMeR == 1);
        model(IsGood,:) = MapToUse(MSClassR(IsGood),:);
        
        fst = mean(model(IsGood,:).*NormDimL2(gmallTPR(IsGood,:),2),2);
        fsl = mean(model(IsGood,:).*NormDimL2(gmallLPR(IsGood,:),2),2);
        
        res(run,n) = mean(fst);
        lres(run,n) = mean(fsl);
        figure(5);
        subplot(121);
        plot(lres(1:run,:)','-','Color',[0.7 0.7 0.7]);
        axis([(out.OptStart - 1) (out.OptEnd + 1) 0 1]);
        title('Learning set');
        ylabel('Fit');
        grid on;
        subplot(122);
        plot(res(1:run,:)','-','Color',[0.7 0.7 0.7]);
        axis([(out.OptStart - 1) (out.OptEnd + 1) 0 1]);
        title('Test set');
        grid on;
    end
    figure(5);
    subplot(121);
    hold off
    plot(lres(1:run,:)','-','Color',[0.7 0.7 0.7]);
    hold on
    ylabel('Fit');
    MeanFit = mean(lres(1:run,:),1);
    plot(out.OptStart:out.OptEnd,MeanFit(1,out.OptStart:out.OptEnd),'-ok','LineWidth',2,'MarkerFaceColor',[0 0 0]);
    axis([(out.OptStart - 1) (out.OptEnd + 1) 0 1]);
    title('Learning set');
    grid on;
    subplot(122);
    hold off
    plot(res(1:run,:)','-','Color',[0.7 0.7 0.7]);
    hold on
    MeanFit = mean(res(1:run,:),1);
    plot(out.OptStart:out.OptEnd,MeanFit(1,out.OptStart:out.OptEnd),'-ok','LineWidth',2,'MarkerFaceColor',[0 0 0]);
    axis([(out.OptStart - 1) (out.OptEnd + 1) 0 1]);
    title('Test set');
    grid on;
end

dummy = inputdlg('Where in the x-axis would you say the plateau of explained variance begins in the test set?','Choose optimal number of microstates');

dummy = str2double(dummy);

if isnan(dummy)
    return
end

out.XValResults = res';
out.nStates = dummy;
ret = ComputeMSMaps(out,[],true);

%[filename, pathname] = uiputfile('*.mat', 'Save Crossvalidation result as');
%if filename == 0
%    return;
%else
%    save(fullfile(pathname,filename),'res');
%end

% --------------------------------------------------------------------
function MenuAnalysis_GFP_Callback(hObject, eventdata, handles)
% hObject    handle to MenuAnalysis_GFP (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
DoTheTanova(hObject, eventdata, handles,true);


% --------------------------------------------------------------------
function View_GFP_Callback(hObject, eventdata, handles)
% hObject    handle to View_GFP (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

out = get(handles.output,'UserData');
out.DoGFP = 1;
set(handles.output,'UserData',out);


if ~isempty(out.GFPPTanova)
    cld();

    Randomizer_ShowTanovaResults(out,handles.figure1);
else
    switch questdlg('Run the analysis now?','Results not yet computed','Yes','No','Yes');
        case 'Yes'
            DoTheTanova(hObject, eventdata, handles,true);
    end
end


% --------------------------------------------------------------------
function MS_SaveMaps_Callback(hObject, eventdata, handles)
% hObject    handle to MS_SaveMaps (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[filename, pathname] = uiputfile('*.asc', 'Save as textfile');
if isequal(filename,0) || isequal(pathname,0)
    return
end
out = get(handles.output,'UserData');

[fp,err] = fopen(fullfile(pathname, filename),'wt');

if fp == -1
    errordlg(err);
    return
end

for i = 1:size(out.MSMaps,1)
    for j = 1:size(out.MSMaps,2)
        if j == size(out.MSMaps,2)
            fprintf(fp,'%f\n',out.MSMaps(i,j));
        else
            fprintf(fp,'%f\t',out.MSMaps(i,j));
        end
    end
end

    
fclose(fp);


% --------------------------------------------------------------------
function T_Mapper_Callback(hObject, eventdata, handles)
% hObject    handle to T_Mapper (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

out = get(handles.figure1,'UserData');
h = Ragu_TMapper(hObject,eventdata,handles,out);

if (isempty(h))
    return
end

dat = get(h,'UserData');

out.LoretaSysDir = dat.LoretaSysDir;
out.LoretaSPINV  = dat.LoretaSPINV;

close(h);
set(handles.figure1,'UserData',out);


% --------------------------------------------------------------------
function Data_LorOptions_Callback(hObject, eventdata, handles)
% hObject    handle to Data_LorOptions (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

dat = get(handles.figure1,'UserData');
h = Ragu_LoretaOptions(hObject, eventdata, handles,dat);

if (isempty(h))
    return
end

out = get(h,'UserData');

dat.LoretaSysDir = out.LoretaSysDir;
dat.LoretaSPINV  = out.LoretaSPINV;

close(h);
cld();
set(handles.figure1,'UserData',dat);


% --------------------------------------------------------------------
function Data_Export_Callback(hObject, eventdata, handles)
% hObject    handle to Data_Export (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

dat = get(handles.figure1,'UserData');

if ~isempty(dat)
    RaguDumpToVision(dat);
end


% --------------------------------------------------------------------
function Options_ResetWindow_Callback(hObject, eventdata, handles)
% hObject    handle to Options_ResetWindow (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

dat = get(handles.figure1,'UserData');
dat.StartFrame = 1;
dat.EndFrame = size(dat.V,4);
set(handles.figure1,'UserData',dat);


% --------------------------------------------------------------------
function FileSave_Callback(hObject, eventdata, handles)
% hObject    handle to FileSave (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

rd = get(handles.output,'UserData');

if isfield(rd,'Modified')
    rd = rmfield(rd,'Modified');
end

save(rd.FileName,'rd');

%cd(pathname);

set(handles.output,'UserData',rd);
set(handles.figure1,'Name',['Ragu: ' rd.FileName]);


% --------------------------------------------------------------------
function MenuFileExpGFP_Callback(hObject, eventdata, handles)
% hObject    handle to MenuFileExpGFP (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

MenuFileExpTanova_Callback(hObject, eventdata, handles,1)


% --------------------------------------------------------------------
function View_DataCursor_Callback(hObject, eventdata, handles)
% hObject    handle to View_DataCursor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

h = datacursormode(handles.figure1);
if strcmp(get(h,'Enable'),'off')
    set(h,'Enable','on');
    set(hObject,'Checked','on');
else
    set(h,'Enable','off');
    set(hObject,'Checked','off');
end


% --------------------------------------------------------------------
function uitoggletool3_OffCallback(hObject, eventdata, handles)
% hObject    handle to uitoggletool3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.View_DataCursor,'Checked','off');

% --------------------------------------------------------------------
function uitoggletool3_OnCallback(hObject, eventdata, handles)
% hObject    handle to uitoggletool3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.View_DataCursor,'Checked','on');


% --------------------------------------------------------------------
function Edit_AddLabel_Callback(hObject, eventdata, handles)
% hObject    handle to Edit_AddLabel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
txt = inputdlg('Enter text');
gtext(txt);

% --------------------------------------------------------------------
function MenuAUCStats_Callback(hObject, eventdata, handles)
% hObject    handle to MenuAUCStats (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

d = get(handles.output,'UserData');

d.TanovaAUC = squeeze(sum(log(d.PTanova(:,:,d.StartFrame:d.EndFrame,:)),3)) ./ (d.EndFrame-d.StartFrame+1) * (-2);
for i = 1:2
    for j = 1:4
        d.PAUCCount(i,j) = squeeze(sum(d.TanovaAUC(i,j,:)>= d.TanovaAUC(i,j,1),3)) / size(d.TanovaAUC,3);
    end
end

set(handles.output,'UserData',d);
Randomizer_ShowHitCountResults(d,2);

% --------------------------------------------------------------------
function MenuAUCStatsGFP_Callback(hObject, eventdata, handles)
% hObject    handle to MenuAUCStatsGFP (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
d = get(handles.output,'UserData');

d.GFPAUC = squeeze(sum(log(d.GFPPTanova(:,:,d.StartFrame:d.EndFrame,:)),3)) ./ (d.EndFrame-d.StartFrame+1) * (-2);
for i = 1:2
    for j = 1:4
        d.PAUCCountGFP(i,j) = squeeze(sum(d.GFPAUC(i,j,:)>= d.GFPAUC(i,j,1),3)) / size(d.GFPAUC,3);
    end
end

set(handles.output,'UserData',d);
Randomizer_ShowHitCountResults(d,4);





% --------------------------------------------------------------------
function Ragu_Help_Import_Callback(hObject, eventdata, handles)
% hObject    handle to Ragu_Help_Import (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Ragu_Import_Help(2);


% --------------------------------------------------------------------
function Ragu_Help_General_Callback(hObject, eventdata, handles)
% hObject    handle to Ragu_Help_General (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Ragu_Import_Help(1);


% --------------------------------------------------------------------
function Ragu_Help_Design_Callback(hObject, eventdata, handles)
% hObject    handle to Ragu_Help_Design (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Ragu_Import_Help(3);


% --------------------------------------------------------------------
function Ragu_Help_Analyze_Callback(hObject, eventdata, handles)
% hObject    handle to Ragu_Help_Analyze (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Ragu_Import_Help(6);


% --------------------------------------------------------------------
function Ragu_Help_Files_Callback(hObject, eventdata, handles)
% hObject    handle to Ragu_Help_Files (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Ragu_Import_Help(4);


% --------------------------------------------------------------------
function Ragu_Help_Options_Callback(hObject, eventdata, handles)
% hObject    handle to Ragu_Help_Options (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Ragu_Import_Help(5);


% --------------------------------------------------------------------
function Ragu_Help_View_Callback(hObject, eventdata, handles)
% hObject    handle to Ragu_Help_View (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Ragu_Import_Help(7);


% --------------------------------------------------------------------
function Tanova_Overall_Callback(hObject, eventdata, handles)
% hObject    handle to Tanova_Overall (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function GFP_Overall_Callback(hObject, eventdata, handles)
% hObject    handle to GFP_Overall (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function MS_Cormat_Callback(hObject, eventdata, handles)
% hObject    handle to MS_Cormat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

out = get(handles.figure1,'UserData');

h = Ragu_Cormat(hObject,eventdata,handles,out);

if (isempty(h))
    return
end

close(h);


% --------------------------------------------------------------------
function Help_Graph_Callback(hObject, eventdata, handles)
% hObject    handle to Help_Graph (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

disp(handles.CurrentView);


% --------------------------------------------------------------------
function Menu_Tools_Callback(hObject, eventdata, handles)
% hObject    handle to Menu_Tools (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Tools_FindOutlier_Callback(hObject, eventdata, handles)
% hObject    handle to Tools_FindOutlier (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

out = get(handles.figure1,'UserData');

h = Ragu_DetectOutlier(hObject,eventdata,handles,out);

if (isempty(h))
    return
end

UseMe = get(h,'UserData');

close(h);

if strcmp(questdlg('Exclude the selected cases?', 'Outlier removal', 'Yes', 'No', 'Yes'),'Yes')
    out.IndFeature(UseMe == false) = nan;
    set(handles.figure1,'UserData',InitResults(out));
end

% --------------------------------------------------------------------
function Tools_RemoveBL_Callback(hObject, eventdata, handles)
% hObject    handle to Tools_RemoveBL (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


out = get(handles.figure1,'UserData');

if ~isfield(out,'OrgData')
    out.OrgData = out.V;
end

Output = Randomizer_Tanova_Options(hObject, eventdata, handles,out,true,true);
out = get(Output,'UserData');
if isempty(out)
    return;
end
close(Output);

out = InitResults(out);
out.V = out.V - repmat(mean(out.V(:,:,:,out.StartFrame:out.EndFrame),4),[1,1,1,size(out.V,4)]);
out.StartFrame = 1;
out.EndFrame   = size(out.V,4);
out.BaselineFramesRemoved = [out.StartFrame out.EndFrame];

set(handles.figure1,'UserData',InitResults(out));
UpdateGUI(handles);

cld();
Randomizer_ShowDataSet2();


% --------------------------------------------------------------------
function ToolsFilter_Callback(hObject, eventdata, handles)
% hObject    handle to ToolsFilter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

out = get(handles.figure1,'UserData');

if ~isfield(out,'OrgData')
    out.OrgData = out.V;
end

if isfield(out,'FilterPars')
    if strcmp(questdlg('Data has previously been filtered!', 'Warning','Cancel', 'Continue', 'Cancel'),'Cancel') == 1
        return
    end
end

DoLoop = true;

fs = 1000 / out.DeltaX;

while DoLoop == true

    answer = inputdlg({'Low cut (Hz), or empty','High cut (Hz), or empty','Order'},'IIR filter parameters');

    if isempty(answer)
        return
    end
 
    out.FilterPars = str2double(answer);
    if isnan(out.FilterPars(3))
        uiwait(msgbox('No valid filter specified','Warning','modal'));
        continue;
    end
    % Low cut there?
    if ~isnan(out.FilterPars(1))
        if ~isnan(out.FilterPars(2))
            % We have a bandpass
            iirpars = fdesign.bandpass('N,F3dB1,F3dB2',out.FilterPars(3),out.FilterPars(1)/fs*2,out.FilterPars(2)/fs*2);
        else
            % We have a low cut
            iirpars = fdesign.highpass('N,F3dB',out.FilterPars(3),out.FilterPars(1)/fs*2);
        end
    else
       % We have no low cut
        if ~isnan(out.FilterPars(2))
            % We have a high cut
            iirpars = fdesign.lowpass('N,F3dB',out.FilterPars(3),out.FilterPars(2)/fs*2);
        else
            % We have crap
            iirpars = nan;
            uiwait(msgbox('No valid filter specified','Warning','modal'));
            continue;
            
        end
    end

    if ~isnumeric(iirpars)
        myfilt = design(iirpars,'IIR');
        [res,w] = freqz(myfilt,1024,fs);
        figure(200);
        plot(w,abs(res));
        switch questdlg('Is this filter ok?', 'Filter response','OK', 'Adjust', 'Cancel','Adjust');
            case 'OK'
                break;
            case 'Cancel'
                close(200);
                return
            otherwise
                DoLoop = true;
        end
    end
end
close(200);


dims = size(out.V);

cld();

h = waitbar(0,'Filtering, please wait...');
            
for s = 1:dims(1)
    for c = 1:dims(2)
        waitbar((s*dims(2) + c)/prod(dims(1:2)),h)
        for i = 1:dims(3);
            % out.V(s,c,i,:) = filtfilt(b,1,squeeze(out.V(s,c,i,:)));
            out.V(s,c,i,:) = filter(myfilt,squeeze(out.V(s,c,i,:)));

        end
    end
end

close(h);
out = InitResults(out);

UpdateGUI(handles);            
set(handles.figure1,'UserData',out);
Randomizer_ShowDataSet2();


% --------------------------------------------------------------------
function ToolsUndo_Callback(hObject, eventdata, handles)
% hObject    handle to ToolsUndo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

out = get(handles.figure1,'UserData');
out.V = out.OrgData;
out = rmfield(out,'OrgData');
if isfield(out,'FilterPars');
    out = rmfield(out,'FilterPars');
end
if isfield(out,'BaselineFramesRemoved')
    out = rmfield(out,'BaselineFramesRemoved');
end

set(handles.figure1,'UserData',InitResults(out));
UpdateGUI(handles);
cld();
Randomizer_ShowDataSet2();


% --------------------------------------------------------------------
function ToolsFilterSpecs_Callback(hObject, eventdata, handles)
% hObject    handle to ToolsFilterSpecs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

out = get(handles.figure1,'UserData');
fs = 1000 / out.DeltaX;

if ~isfield(out,'FilterPars');
    out.FilterPars = nan(3,1);
end

if ~isnan(out.FilterPars(1))
    msg{1} = sprintf('Low cut: %f', out.FilterPars(1));
else
    msg{1} = 'No low cut applied.';
end
    
if ~isnan(out.FilterPars(2))
    msg{2} = sprintf('High cut: %f', out.FilterPars(2));
else
    msg{2} = 'No high cut applied.';
end

if ~isnan(out.FilterPars(3))
    msg{3} = sprintf('Filter order: %i', out.FilterPars(3));
else
    msg{3} = 'No filter order defined.';
end

if isfield(out,'BaselineFramesRemoved')
    msg{4} = 'Baseline correction has been applied';
else
    msg{4} = 'No Baseline correction has been applied';
end

uiwait(msgbox(msg,'Filter settings','modal'));
 
