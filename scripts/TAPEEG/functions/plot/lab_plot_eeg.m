% script to view eeg/meg data, using function 'eegplot' of EEGLAB
%
% lab_plot_eeg(data,header)
%
% written by F. Hatz 2013

function [data,header,cfg] = lab_plot_eeg(data,header,cfg,dostore)

storeeeg = false;

if ~exist('data','var') & ~exist('header','var')
    [data,header,cfg] = lab_read_data;
    if isempty(header)
        return
    end
end

if ~exist('cfg','var')
    cfg = [];
end
if ~exist('header','var')
    header = lab_create_header(data);
end

if ~isfield(header,'EEG_file')
    header.EEG_file = 'EEG_file';
end
for i = 1:size(header.channels,1)
    eloc(1,i).labels = header.channels(i,:); %#ok<AGROW>
    if isfield(header,'badchans') & find(header.badchans == i)
        eloc(1,i).badchan = 1; %#ok<AGROW>
    else
        eloc(1,i).badchan = 0; %#ok<AGROW>
    end
    eloc(1,i).index = i; %#ok<AGROW>
end
if isfield(header,'events') & isfield(header.events,'POS')
   events = [];
   for i = 1:size(header.events.POS,2)
       events(1,end+1).type = header.events.TYP{1,i}; %#ok<AGROW>
       %events(i).condition = header.events.TYP{1,i};
       events(1,end).latency = double(header.events.POS(1,i)-1);
       events(1,end).duration = double(header.events.DUR(1,i));
       events(1,end).index=length(events);
       events(1,end).proc='none';
   end
else
    events = [];
end

set(0,'units','pixels');
pos = get(0,'ScreenSize');
pos = [50 80 pos(3)-100 pos(4)-150];
eegplot(data,'srate',header.samplingrate,'dispchans',20,'winlength',10, ...
    'xgrid','on','title',header.EEG_file,'eloc_file',eloc,'events',events, ...
    'ploteventdur','on','butlabel','Store','position',pos, ...
    'ctrlselectcommand',{'lab_ctrldowncom('' '',''off'');','eegplot(''defmotioncom'',gcbf);',''}, ...
    'command','assignin(''caller'',''reject'',TMPREJ)');
    % 'ctrlselectcommand',{'lab_ctrldowncom('' '',''off'');','eegplot(''defmotioncom'',gcbf);',''}, ...
uiwait;

try
    reject = evalin('base','TMPREJ');
    dostore = 1;
catch %#ok<CTCH>
    reject = [];
end
try
    eventsout = evalin('base','EVENTS');
catch %#ok<CTCH>
    eventsout = [];
end
try
    badchans = evalin('base','BAD');
catch %#ok<CTCH>
    badchans = [];
end
if ~isempty(eventsout)
    events = [];
    for i = 1:length(eventsout)
        events.POS(1,i) = int64(eventsout(i).latency+1);
        events.DUR(1,i) = int64(eventsout(i).duration);
        events.OFF(1,i) = int64(0);
        events.TYP{1,i} = eventsout(i).type;
    end
    header.events = events;
    storeeeg = true;
end
if ~isempty(reject)
    if ~isfield(header,'events') | ~isfield(header.events,'POS')
        header.events.POS = [];
        header.events.DUR = [];
        header.events.OFF = [];
        header.events.TYP = {};
    end
    for i = 1: size(reject,1)
        Estart = round(reject(i,1));
        Edur = (round(reject(i,2)) - Estart) + 1;
        if Estart < 1
            Estart = 1;
        end
        if Edur < 1
            Edur = 1;
        end
        header.events.POS(1,end+1) = int64(Estart);
        header.events.DUR(1,end+1) = int64(Edur);
        header.events.OFF(1,end+1) = int64(0);
        header.events.TYP{1,end+1} = 'BAD';
    end
    [~,Esort] = sort(header.events.POS);
    header.events.POS = header.events.POS(1,Esort);
    header.events.DUR = header.events.DUR(1,Esort);
    header.events.OFF = header.events.OFF(1,Esort);
    header.events.TYP = header.events.TYP(1,Esort);
    
    storeeeg = true;
end
if ~isempty(badchans)
    for i = 1:length(badchans)
        bad(i) = badchans(i).badchan; %#ok<AGROW>
    end
    clearvars badchans
    badchans = find(bad==1);
    if ~isempty(badchans)
        header.badchans = badchans;
        header.goodchans = setdiff(header.numdatachannels,badchans);
        header.goodchans = header.goodchans(:)';
    end
end
if ~exist('data','var')
    data = lab_read_data(fullfile(cfg.EEG_filepath,cfg.EEG_file));
end

if storeeeg == true & exist('cfg','var') & isfield(cfg,'EEG_file')
    if exist('dostore','var') & dostore == 1
        answer = 'Yes';
    else
        answer = questdlg('Store modifications?','Save data','No','Yes','Yes');
    end
    if strcmp(answer,'Yes')
        lab_save_data(data,header,fullfile(cfg.EEG_filepath,cfg.EEG_file));
        [~,~,~,FilenameS] = lab_filename(cfg.EEG_file);
        if exist(fullfile(cfg.EEG_filepath,[FilenameS '_exclude~.txt']),'file')
            delete(fullfile(cfg.EEG_filepath,[FilenameS '_exclude~.txt']));
        end
        if exist(fullfile(cfg.EEG_filepath,[FilenameS '_exclude.txt']),'file')
            delete(fullfile(cfg.EEG_filepath,[FilenameS '_exclude.txt']));
        end
    end
end

disp('')
if nargout == 0
    data= [];
    header = [];
    cfg = [];
end

end