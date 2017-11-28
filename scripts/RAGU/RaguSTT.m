function varargout = RaguSTT(varargin)
% RAGUSTT M-file for RaguSTT.fig
%      RAGUSTT, by itself, creates a new RAGUSTT or raises the existing
%      singleton*.
%
%      H = RAGUSTT returns the handle to a new RAGUSTT or the handle to
%      the existing singleton*.
%
%      RAGUSTT('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in RAGUSTT.M with the given input arguments.
%
%      RAGUSTT('Property','Value',...) creates a new RAGUSTT or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before RaguSTT_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to RaguSTT_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Copyright 2009-2011 Thomas Koenig
% distributed under the terms of the GNU AFFERO General Public License

% Edit the above text to modify the response to help RaguSTT

% Last Modified by GUIDE v2.5 09-Sep-2010 14:46:13

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @RaguSTT_OpeningFcn, ...
                   'gui_OutputFcn',  @RaguSTT_OutputFcn, ...
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


% --- Executes just before RaguSTT is made visible.
function RaguSTT_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to RaguSTT (see VARARGIN)

% Choose default command line output for RaguSTT
in = varargin{4};
handles.output = hObject;

rn{1} = 'Main effect';
rn{2} = in.strF1;
rn{3} = in.strF2;
rn{4} = [in.strF1 ' x ' in.strF2];

set(handles.uitable1,'RowName',rn');

cn{1} = 'Main effect';
cn{2} = in.IndName;
set(handles.uitable1,'ColumnName',cn');
set(handles.uitable1,'Data',in.stt');

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes RaguSTT wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = RaguSTT_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
