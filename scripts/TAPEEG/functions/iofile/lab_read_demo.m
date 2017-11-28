% Demo read file as template for own import scripts
%
% [data,header,cfg] = lab_read_demo(filename)
%
% written by F. Hatz

function [data,header,cfg] = lab_read_demo(filename)

if ~exist('filename','var')
    [filepath,filenameS] = uigetfile('*.*','select EEG-file');
    filename = fullfile(filepath,filenameS);
else
    [filenameS,filepath] = lab_filename(filename);
end

% Read data
data = [];
%    some function to read the data
%    data must be formatted channels x timeframes
%    auxillary channels should be last channels

% Read header information
samplingrate =  []; % number
auxchans =  []; % number of auxillary channels (last channels)
labels =  []; % textcell-array (1 x channels)
datestr =  []; % string with date&time 'yyyy-mm-dd hh:mm:ss.miliseconds'
PatName =  []; % string with patient name
PatDOB =  []; % string with patient day-of-birth 'yyyy-mm-dd'
PatSex =  []; % string with patient gender 'F' or 'M'
PatID =  []; % string with patient ID

header.version=0;
header.samplingrate = samplingrate;
header.numchannels = size(data,1);
header.numdatachannels=size(data,1) - auxchans;
header.numauxchannels = hdr.auxchans;
header.numtimeframes=size(data,2);
header.channels = char(labels(1:header.numchannels,1)');
[header.year,header.month,header.day, header.hour, header.minute, header.second]= datevec(datestr);
header.millisecond = (header.second - floor(header.second))*1000;
header.second = floor(header.second);
header.subject.ID = PatID;
header.subject.name = PatName;
[header.subject.year,header.subject.month,header.subject.day] = datevec(PatDOB,'yyyy-mm-dd');
header.subject.sex = PatSex;

% Read events
eventsPOS =  []; % vector (1 x number events) with position in timeframes
eventsDUR =  []; % vector (1 x number events) with duration in timeframes
eventsTYP =  []; % text-cell-array (1 x number events) with names

header.events.POS = int64(eventsPOS);
header.events.DUR = int64(eventsDUR);
header.events.OFF = int64(zeros(1,size(eventsPOS,2)));
header.events.TYP = eventsTYP;

% read electrodes locations
LOCS =  []; % filepath to electrodes-file
header.locs = lab_read_locs(LOCS);

header.datatype =  []; % 'eeg' or 'meg'
header.ref_chan =  []; % string or number(s)
header.EEG_file = filenameS;
header.EEG_filepath = filepath;

cfg.EEG_file = filenameS;
cfg.EEG_filepath = filepath;

return
