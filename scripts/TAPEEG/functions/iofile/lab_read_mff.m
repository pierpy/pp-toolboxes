% Function to read Netstation mff
%
% [data,header,cfg] = lab_read_mff(filename,cfg)
%
% written by F. Hatz / H. Bousleiman 2013

function [data,header,cfg] = lab_read_mff(Filename,cfg,nodata,segment)

if ~exist('segment','var')
    segment = [];
end
if ~exist('nodata','var')
    nodata = false;
end
if ~exist('cfg','var')
    cfg = [];
end
if ~exist('Filename','var')
    [~,Filename] = uigetfile('*.bin','select signal1.bin in mff-folder');
end

% Convert Filename to filepath
tmp = strfind(Filename,filesep);
if strcmp(Filename(end-3:end),'.mff') & exist(fullfile(Filename(1:end-4),'signal1.bin'),'file')
    Filepath = Filename(1:end-4);
elseif exist(fullfile(Filename(1:tmp(end)),'signal1.bin'),'file')
    Filepath = Filename(1:tmp(end)-1);
elseif exist(fullfile(Filename,'signal1.bin'),'file')
    Filepath = Filename;
else
    data = [];
    header = [];
    disp('Abort, no valid mmf-file selected')
    return
end
clearvars tmp

% delete some weird mac stuff
delete(fullfile(Filepath,'._*'));

% Read data
warning off %#ok<WNOFF>
[hdr] = read_mff_header(Filepath);
if strcmp(hdr.orig.epochType,'cnt')
    format = 'sample';
elseif strcmp(hdr.orig.epochType,'seg')
    format = 'seg';
else
    data = [];
    header = [];
    disp('Error reading mff (unknown filetype)')
    return
end

% correct for old mff versions
InfoObj = mff_getObject(com.egi.services.mff.api.MFFResourceType.kMFF_RT_Info,'info.xml',Filepath);
if InfoObj.getMFFVersion == 0
    hdr.nSamples = hdr.nSamples / hdr.Fs;
    if strcmp(format,'sample')
        hdr.orig.epochBeginSamps = hdr.orig.epochBeginSamps / hdr.Fs;
        hdr.orig.epochNumSamps = hdr.orig.epochNumSamps / hdr.Fs;
        hdr.orig.epochTime0 = hdr.orig.epochTime0 / hdr.Fs;
    end
end

try
    if nodata == false
        if ~isempty(segment) & length(segment) > 1
            if segment(1) > hdr.nSamples-1
                segment(1) = hdr.nSamples-1;
            end
            if segment(2) > hdr.nSamples
                segment(2) = hdr.nSamples;
            end
            [data] = read_mff_data(Filepath,'sample',segment(1),segment(2),1:hdr.nChans,hdr);
        else
            [data] = read_mff_data(Filepath,'sample',1,hdr.nSamples,1:hdr.nChans,hdr);
        end
        data = double(data); % TAPEEG needs double precision!
    else
        data = zeros(hdr.nChans,1);
    end
catch %#ok<CTCH>
    data = [];
    header = [];
    disp(['Error reading mff (' regexprep(Filepath,filesep,'//') ')'])
    return
end

% define lag
[lag,cfg] = lab_get_egi_lag(hdr.Fs,Filepath,cfg);
if lag > 0
    disp(['    correct external events and analog channels with fixed lag of ' num2str(lag)])
end

% Read events
[eventtmp] = read_mff_event(Filepath,hdr);
if ~isempty(eventtmp)
    eventtmp = struct2cell(eventtmp');
    events.TYP = eventtmp(1,:);
    events.POS = int64(cell2mat(eventtmp(2,:)));
    events.DUR = int64(cell2mat(eventtmp(5,:)));
    events.OFF = int64(cell2mat(eventtmp(4,:)));
    if size(eventtmp,1) >= 9
        external = eventtmp(9,:);
    end
    if isempty(events.OFF)
        events.OFF = int64(zeros(1,length(events.POS)));
    end
    header.events = events;
    
    % convert epoch/segment info to markers for TAPEEG
    eventstmp = events;
    tmp = find(strcmp(events.TYP,'break cnt'));
    tmp2 = setdiff(1:size(events.TYP,2),tmp);
    events.TYP = events.TYP(1,tmp2);
    events.POS = events.POS(1,tmp2);
    events.DUR = events.DUR(1,tmp2);
    events.OFF = events.OFF(1,tmp2);
    if exist('external','var') & length(external) >= max(tmp2)
        external = external(tmp2);
        for i = 1:length(events.POS)
            if external{i} == 1
                events.POS(i) = events.POS(i) + int64(lag);
            end
        end
    end
    clearvars tmp2
    for i = tmp
        if eventstmp.POS(1,i) > 2 * hdr.Fs | eventstmp.DUR(1,i) > 6 * hdr.Fs
            events.TYP = [events.TYP {'SegStart','SegStop'}];
            events.POS = [events.POS eventstmp.POS(1,i) (eventstmp.POS(1,i) + eventstmp.DUR(1,i) - 1)];
            events.DUR = [events.DUR 1 1];
            events.OFF = [events.OFF 0 0];
        end
    end
    [events.POS,tmp] = sort(events.POS);
    events.TYP = events.TYP(tmp);
    events.DUR = events.DUR(tmp);
    events.OFF = events.OFF(tmp);
    clearvars tmp
    header.events = events;
end
clearvars eventtmp;

% Set auxillary channels
aux = ones(1,size(data,1),1);
for i = 1:size(data,1)
    if strcmp(hdr.chantype{i},'eeg')
        aux(i) = 0;
    end
end
Iaux = find(aux == 1);
if ~isempty(Iaux)
    hdr.auxchans = length(Iaux);
    Iaux = Iaux(:)';
    Ieeg = setdiff(1:size(data,1),Iaux);
    if min(Iaux) < max(Ieeg)
        data = data([Ieeg Iaux],:);
        hdr.label = hdr.label(1,[Ieeg Iaux]);
        hdr.chantype = hdr.chantype(1,[Ieeg Iaux]);
        hdr.chanunit = hdr.chanunit(1,[Ieeg Iaux]);
    end
    tmp = find(strcmpi(hdr.chantype,'ECG'));
    if ~isempty(tmp)
        hdr.label{tmp(1)} = 'ECG';
    end
    tmp = find(strcmpi(hdr.chantype,'EOG'));
    if ~isempty(tmp)
        hdr.label{tmp(1)} = 'EOG';
    end
    if lag > 0 & size(data,2) > 1
        data(Iaux,:) = [repmat(data(Iaux,1),1,lag) data(Iaux,1:end-lag)];
    end
else
    hdr.auxchans = 0;
end

clearvars Iaux Ieeg aux

% Collect header information
header.samplingrate = hdr.Fs;
header.numchannels = size(data,1);
header.version = 0;
header.numdatachannels=size(data,1) - hdr.auxchans;
header.numauxchannels = hdr.auxchans;
if size(data,2) == 1
    header.numtimeframes = hdr.nSamples;
else
    header.numtimeframes = size(data,2);
end
header.channels = char(hdr.label(1,1:header.numchannels)');
hdr = rmfield(hdr,'orig');
header.orighdr = hdr;

% read record time from info.xml
recordTime = mff_getObject(com.egi.services.mff.api.MFFResourceType.kMFF_RT_Info,'info.xml',Filepath);
recordTime = char(recordTime.getRecordTime);
tmp = strfind(recordTime,'T');
recordTime2 = recordTime(tmp+1:end);
recordTime = recordTime(1:tmp-1);
tmp = strfind(recordTime2,'+');
if ~isempty(tmp)
    recordTime2 = recordTime2(1:tmp-1);
end
tmp = strfind(recordTime2,'-');
if ~isempty(tmp)
    recordTime2 = recordTime2(1:tmp-1);
end
recordTime = [recordTime ' ' recordTime2];
[header.year,header.month,header.day, header.hour, header.minute, header.second]= datevec(recordTime);
header.millisecond = (header.second - floor(header.second))*1000;
header.second = floor(header.second);
clearvars recordTime recordTime2 tmp

% Read subject info from subject.xml
try %#ok<TRYNC>
    subject = xmltools(fullfile(Filepath,'subject.xml'));
    for i = 1:length(subject.child(1,2).child(1,1).child)
        value1 = subject.child(1,2).child(1,1).child(1,i).child(1,1).value;
        value2 = subject.child(1,2).child(1,1).child(1,i).child(1,2).value;
        value2 = regexprep(value2, '[<>]', '');
        switch value1
            case 'Patient-ID'
                header.subject.ID = value2;
            case 'localIdentifier'
                header.subject.ID = value2;
            case 'Nachname'
                header.subject.name = value2;
            case 'Vorname'
                if isfield(header.subject,'name')
                    header.subject.name = [value2 '_' header.subject.name];
                end
            case 'Geburtsdatum'
                [header.subject.year,header.subject.month,header.subject.day] = datevec(value2,'yyyy-mm-dd');
            case 'Date of Birth'
                [header.subject.year,header.subject.month,header.subject.day] = datevec(value2,'yyyy-mm-dd');
            case 'Geschlecht'
                if strcmp(value2,'false')
                    header.subject.sex = 'M';
                else
                    header.subject.sex = 'F';
                end
            case 'Gender'
                if strcmp(value2,'false')
                    header.subject.sex = 'M';
                else
                    header.subject.sex = 'F';
                end
            case 'Dominante Hand'
                header.subject.hand = value2;
        end
    end
end
warning on %#ok<WNON>

% Take subject/time info from Filename (when missing)
tmp = strfind(Filepath,filesep);
if length(tmp) > 1
    if tmp(end) == length(Filepath)
        patinfo = textscan(Filepath(tmp(end-1)+1:end-1),'%s');
    else
        patinfo = textscan(Filepath(tmp(end)+1:end),'%s');
    end
    patinfo = patinfo{1,1};
    if size(patinfo,1) > 1
        for i = 1:size(patinfo,1)
            if isnan(str2double(patinfo{i,1}))
                if ~isfield(header,'subject') || ~isfield(header.subject, 'name') || length(header.subject.name) < 3
                    header.subject.name =  regexprep(patinfo{i,1},',','');
                    flagname = 1; %#ok<NASGU>
                elseif i == 2 & exist('flagname','var')
                    header.subject.name = [regexprep(patinfo{i,1},',','') '_' header.subject.name];
                end
            else
                if length(patinfo{i,1}) == 8 && header.year == 0
                    header.year = str2num(patinfo{i,1}(1:4)); %#ok<ST2NM>
                    header.month = str2num(patinfo{i,1}(5:6)); %#ok<ST2NM>
                    header.day = str2num(patinfo{i,1}(7:8)); %#ok<ST2NM>
                end
                if length(patinfo{i,1}) == 4 && header.hour == 0
                    header.hour = str2num(patinfo{i,1}(1:2)); %#ok<ST2NM>
                    header.minute = str2num(patinfo{i,1}(3:4)); %#ok<ST2NM>
                end
            end
        end
    end
end
clearvars patinfo

if isfield(hdr,'locs')
    header.locs = hdr.locs;
    if sum(header.locs.x)^2 > sum(header.locs.y)^2
        tmp = header.locs.x;
        header.locs.x = -header.locs.y;
        header.locs.y = tmp;
        clearvars tmp
    end
    [th,radius] = cart2pol(header.locs.x,header.locs.y,header.locs.z);
    header.locs.radius = radius;
    header.locs.theta = th;
    if header.numchannels > header.numdatachannels & header.numauxchannels > 0
        header.locs.aux = header.numauxchannels;
    end
end

% check for corrupt data
tmp = find(max(isnan(data),[],2));
if ~isempty(tmp)
    if ~isfield(cfg,'EXTRA') | ~isfield(cfg.EXTRA,'removeNaN')
        if ~exist('cfg','var') | ~isfield(cfg,'MAIN') | ~isfield(cfg.MAIN,'auto') | cfg.MAIN.auto == 0
            cfg.EXTRA.removeNaN = true;
            Prompt(1,1:2) = {['Channels ' num2str(tmp) ' are invalid'],''};
            Formats(1).type = 'text';
            Formats(2,1).type = 'none';
            Prompt(2,1:2)={'Set invalid channels to zero','removeNaN'};
            Formats(3,1).type = 'check';
            cfg.EXTRA = inputsdlg(Prompt,'Correct invalid channels',Formats,cfg.EXTRA);
        else
            cfg.EXTRA.removeNaN = true;
        end
    end
    if cfg.EXTRA.removeNaN == true
        data(tmp,:) = 0;
        tmp = tmp(tmp<=header.numdatachannels);
        header.badchans = tmp;
        header.goodchans = setdiff((1:header.numdatachannels),header.badchans);
        header.goodchans = header.goodchans(:)';
    end
end
