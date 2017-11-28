% Helper function for lab_read_data
%
% written by F. Hatz 2012

function [data,header,cfg] = lab_read_data_ext(data,header,cfg,novrb)

if ~exist('novrb','var')
    novrb = false;
end
    
% Get filename and path
Filename = cfg.EEG_file;
Filepath = cfg.EEG_filepath;
[~,~,FormatS,FilenameS] = lab_filename(Filename);

% Correct for missing numdatachannels and numauxchannels
if ~isfield(header,'numchannels')
    header.numchannels = size(data,1);
end
if ~isfield(header,'numdatachannels') & ~isfield(header,'numauxchannels')
    header.numdatachannels = header.numchannels;
    header.numauxchannels = 0;
elseif ~isfield(header,'numdatachannels')
    header.numdatachannels=header.numchannels-header.numauxchannels;
end

% Read locs if available
if exist(fullfile(Filepath,[FilenameS '.els']),'file')
    LOCS = lab_read_locs(fullfile(Filepath,[FilenameS '.els']));
elseif exist(fullfile(Filepath,[FilenameS '.elc']),'file')
    LOCS = lab_read_locs(fullfile(Filepath,[FilenameS '.elc']));
elseif exist(fullfile(Filepath,[FilenameS '.xyz']),'file')
    LOCS = lab_read_locs(fullfile(Filepath,[FilenameS '.xyz']));
elseif exist(fullfile(Filepath,[FilenameS '.sfp']),'file')
    LOCS = lab_read_locs(fullfile(Filepath,[FilenameS '.sfp']));
elseif exist(fullfile(Filepath,[FilenameS '.spi']),'file')
    LOCS = lab_read_locs(fullfile(Filepath,[FilenameS '.spi']));
elseif exist(fullfile(Filepath,'electrodes.els'),'file')
    LOCS = lab_read_locs(fullfile(Filepath,'electrodes.els'));
elseif exist(fullfile(Filepath,'electrodes.elc'),'file')
    LOCS = lab_read_locs(fullfile(Filepath,'electrodes.elc'));
elseif exist(fullfile(Filepath,'electrodes.sfp'),'file')
    LOCS = lab_read_locs(fullfile(Filepath,'electrodes.sfp'));
elseif exist(fullfile(Filepath,'electrodes.xyz'),'file')
    LOCS = lab_read_locs(fullfile(Filepath,'electrodes.xyz'));
end

if ~isfield(header,'channels')
    if exist('LOCS','var') & isfield(LOCS,'aux') & size(data,1) == (size(LOCS.x,2) + LOCS.aux)
        labels = LOCS.labels;
        for i = 1:LOCS.aux
            labels(1,end+1) = cellstr(['AUX' num2str(i)]); %#ok<AGROW>
        end
        header.channels = char(labels');
        clearvars labels
    else
        labels = cell(header.numchannels,1);
        for i = 1:header.numchannels
           labels{i,1} = ['E' num2str(i)];
        end
        header.channels = char(labels);
    end
end
if isfield(cfg,'LOC_file')
    cfg = rmfield(cfg,'LOC_file');
end

% Look for ECG and AUX channels
MaxChan = [];
header.eog_ch = [];
if size(header.channels,2) > 2
    for i = 1:size(header.channels,1)
        if strcmp(header.channels(i,1:3),'ECG') | strcmp(header.channels(i,1:3),'EKG')
            header.ecg_ch = i;
            MaxChan = min(MaxChan,i);
        end
        if strcmp(header.channels(i,1:3),'EOG')
            header.eog_ch = [header.eog_ch i];
            MaxChan = min(MaxChan,i);
        end
        if strcmp(header.channels(i,1:3),'AUX') | strcmp(header.channels(i,1:3),'MIS') | strcmp(header.channels(i,1:3),'ADC')
            MaxChan = min(MaxChan,i);
        end
        
    end
    if ~isempty(MaxChan) & MaxChan <= header.numdatachannels
        header.numdatachannels = MaxChan - 1;
        header.numauxchannels = header.numchannels - header.numdatachannels;
    end
end
clearvars I MaxChan

% Read EEGinfo -- *.txt included for old eeginfo
if exist(fullfile(Filepath,[FilenameS '.txt']),'file') & ~strcmp(FormatS,'txt')
    try %#ok<TRYNC>
        header = lab_read_eeginfo(fullfile(Filepath,[FilenameS '.txt']),header);
    end
end
if exist(fullfile(Filepath,[FilenameS '.info']),'file')
    try %#ok<TRYNC>
        header = lab_read_eeginfo(fullfile(Filepath,[FilenameS '.info']),header,novrb);
    end
end

% Read individual frequency bands
if exist(fullfile(Filepath,[FilenameS '.ifreq']),'file')
    try %#ok<TRYNC>
        load(fullfile(Filepath,[FilenameS '.ifreq']),'-mat');
        if exist('IFREQ','var')
            header.IFREQ = IFREQ;
        end
    end
end

% Read additional marker files (*.mrk , *.edf)
if ~exist('cfg','var')
    cfg = [];
end
[header,skipprocessing] = lab_read_markers(fullfile(Filepath,Filename),header,novrb);
if skipprocessing == 1
    data = [];
    header = [];
    return
end

% correct for spaces in markers
if isfield(header,'events') & isfield(header.events,'POS')
    for i = 1:length(header.events.TYP)
        header.events.TYP{i} = regexprep(header.events.TYP{i},' ','');
    end
end

% Control number of LOCS
if exist('LOCS','var')
    if size(LOCS.x,2) == header.numdatachannels
        header.locs = LOCS;
    elseif isfield(LOCS,'aux') & LOCS.aux > 0 & size(LOCS.x,2) + LOCS.aux == header.numdatachannels & ...
            header.numchannels == header.numdatachannels
    elseif novrb == false
        disp('     Locs ignored (different number of electrodes)')
    end
    clearvars LOCS
    
    % Read digitizer info if available
    if exist(fullfile(Filepath,[FilenameS '_digits.spi']),'file')
        tmp = lab_read_spi(fullfile(Filepath,[FilenameS '_digits.spi']));
        header.locs.digits(:,1) = tmp.x';
        header.locs.digits(:,2) = tmp.y';
        header.locs.digits(:,3) = tmp.z';
    end
end

% Set locs
if isfield(cfg,'LOCS') & isfield(cfg.LOCS,'matchlocs') & cfg.LOCS.matchlocs == true
    [data,header] = lab_match_locs_eeg(data,header,cfg);
elseif isfield(cfg,'LOCS') & isfield(cfg.LOCS,'locs') & ~isempty(cfg.LOCS.locs)& ...
        (~isfield(cfg.LOCS,'filelocs') | cfg.LOCS.filelocs == false) & ...
        header.numchannels >= length(cfg.LOCS.locs.x)
    if novrb == false
        disp('Replace locs by default locs')
    end
    header.locs = cfg.LOCS.locs;
end

% add extra channels modifications
if isfield(cfg,'EXTRA') & ~isempty(cfg.EXTRA)
    [data,header,cfg] = lab_extrachannels(data,header,cfg);
end

% add reference channel
if isfield(cfg,'ADDREF') & isfield(cfg.ADDREF,'name') & ~isempty(cfg.ADDREF.name)
    [header,data,cfg] = lab_add_refchan(header,data,cfg);
end

% Correct for missing ref_chan
if ~isfield(header,'ref_chan')
    header.ref_chan=[];
end

% Correct for missing ecg_ch
if ~isfield(header,'ecg_ch')
    header.ecg_ch = 0;
end

if exist(fullfile(Filepath,[FilenameS '.mat']),'file') & strcmp(FilenameS(end-2:end),'ICA')
    headerICA = read_ICAheader(fullfile(Filepath,[FilenameS '.mat']));
    if ~isempty(headerICA)
        header.W = headerICA.W;
        header.ICAchans = headerICA.ICAchans;
        header.ICAref_chan = headerICA.ref_chan;
        if isfield(headerICA,'numdatachannels')
            header.ICAmaxchan = length(headerICA.numdatachannels);
        else
            header.ICAmaxchan = length(headerICA.numchannels);
        end
    end
elseif exist(fullfile(Filepath,[FilenameS '.w']),'file')
    header.W = importdata(fullfile(Filepath,[FilenameS '.w']));
end

if exist(fullfile(Filepath,[FilenameS '_exclude.txt']),'file') & strcmp(FilenameS(end-2:end),'ICA')
    fid = fopen(fullfile(Filepath,[FilenameS '_exclude.txt']));
    badactivations = textscan(fid,'%f');
    badactivations = badactivations{1}(:)';
    fclose(fid);
    clearvars fid
    if ~isempty(badactivations)
        disp(['     read bad activations from ' fullfile(Filepath,[FilenameS '_exclude.txt'])])
        header.badchans = badactivations;
        header.goodchans = setdiff(1:header.numdatachannels,header.badchans);
        header.goodchans = header.goodchans(:)';
    end
end

if exist(fullfile(Filepath,[FilenameS '_segment.vrb']),'file')
    header = lab_read_segments(fullfile(Filepath,[FilenameS '.sef']),header);
end

% get microstates
header = get_microstates(data,header);

% Correct for missing datatype
if ~isfield(header,'datatype') | isempty(header.datatype)
    header.datatype = 'eeg';
end

% Store Filename and path in header
header.EEG_file = Filename;
header.EEG_filepath = Filepath;
if ~exist(fullfile(Filepath,[FilenameS '.info']),'file')
    lab_write_eeginfo(fullfile(Filepath,[FilenameS '.info']),header);
end

end

function header = read_ICAheader(Filename)
   disp(['    read ICA-info from ' Filename])
   load(Filename);
   if exist('header','var') & exist('W','var')
       header.W = W;
   else
       header = [];
   end
   if exist('ICAchans','var')
       header.ICAchans = ICAchans;
   else
       header.ICAchans = header.goodchans;
   end
end

function header = get_microstates(data,header)
    Mvalues = [];
    Cvalues = [];
    for i = header.numdatachannels+1:header.numchannels
        if size(header.channels,2) >= 5 & strcmpi(header.channels(i,1:5),'Micro')
            tmp = str2num(header.channels(i,6:end)); %#ok<ST2NM>
            if ~isempty(tmp)
                Mvalues = [Mvalues cat(1,tmp,i)]; %#ok<AGROW>
            end
            clearvars tmp
        elseif size(header.channels,2) >= 4 & strcmpi(header.channels(i,1:4),'Corr')
            tmp = str2num(header.channels(i,5:end)); %#ok<ST2NM>
            if ~isempty(tmp)
                Cvalues = [Cvalues cat(1,tmp,i)]; %#ok<AGROW>
            end
            clearvars tmp
        end
    end
    if isempty(Mvalues)
        return
    end
    if size(Mvalues,2) == 1
        header.microstates = Mvalues(1);
    end
    if size(Cvalues,2) == 1
        disp('     read microstates correlations')
        header.CORR = data(Cvalues(2,1),:);
    end
end