% Find good / bad channels
%
% [header,cfg] = lab_find_good(data,header,cfg)
%
% data     = matrix (chans x timeframes)
% header   = output of lab_read_data
% cfg      = structure with config (optional)
%
% written by F. Hatz 2012

function [header,cfg] = lab_find_good(data,header,cfg)
disp('   Detect good channels')

if ~exist('header','var') | ~isfield(header,'samplingrate')
    header = lab_create_header(data);
end
if ~exist('cfg','var')
    cfg = [];
end
if ~isfield(header,'numdatachannels')
    header.numdatachannels = size(data,1);
end
if ~isfield(cfg,'EEG_file') & isfield(header,'EEG_file')
    cfg.EEG_file = header.EEG_file;
    cfg.EEG_filepath = header.EEG_filepath;
end
if ~isfield(header,'EEG_file')
   header.EEG_file = cfg.EEG_file;
   header.EEG_filepath = cfg.EEG_filepath;
end

if ~isfield(cfg,'BADELEC') | ~isfield(cfg.BADELEC,'percentbad')
    if ~isfield(cfg,'BADELEC')
        cfg.BADELEC = [];
    end
    [cfg,skipprocessing] = lab_set_find_good(cfg,header);
    if skipprocessing == 1
        return
    end
end

[badelectrodes,bad,cfg.BADELEC] = lab_detect_bad(data,header,cfg.BADELEC,cfg);
bad.All = badelectrodes(:)';
header.bad = bad;
clearvars badelectrodes

header.badchans = header.bad.All;
if isfield(cfg,'exclude')
    header.badchans = setdiff(header.badchans,cfg.exclude);
end
if isfield(header,'ref_chan') & isnumeric(header.ref_chan)
    header.badchans = setdiff(header.badchans,header.ref_chan);
end
if isfield(cfg,'exclude')
    header.goodchans = setdiff(1:header.numdatachannels,union(header.bad.All,cfg.exclude));
else
    header.goodchans = setdiff(1:header.numdatachannels,header.bad.All);
end
header.badchans = header.badchans(:)';
header.goodchans = header.goodchans(:)';

if isfield(cfg.BADELEC,'markbad') & cfg.BADELEC.markbad == true & isfield(bad,'epochs');
    if ~isfield(header,'events') | ~isfield(header.events,'POS')
        header.events.POS = [];
        header.events.OFF = [];
        header.events.DUR = [];
        header.events.TYP = {};
    else
        tmp = find(strcmp(header.events.TYP,'BAD'));
        if ~isempty(tmp)
            tmp = setdiff(1:length(header.events.POS),tmp);
            header.events.POS = header.events.POS(tmp);
            header.events.DUR = header.events.DUR(tmp);
            header.events.TYP = header.events.TYP(tmp);
            header.events.OFF = header.events.OFF(tmp);
        end
    end
    duration = ceil(cfg.BADELEC.length * header.samplingrate);
    if ~isfield(cfg.BADELEC,'markbadValue')
        markbadValue = 0.25;
    else
        markbadValue = cfg.BADELEC.markbadValue /100;
    end
    bad = mean(bad.epochs,1);
    bad = (find(bad>markbadValue)-1) * duration + 1;
    
    header.events.POS = [header.events.POS int64(bad)];
    header.events.DUR = [header.events.DUR repmat(duration,1,length(bad))];
    header.events.OFF = [header.events.OFF zeros(1,length(bad))];
    header.events.TYP = [header.events.TYP repmat(cellstr('BAD'),1,length(bad))];
    [~,tmp] = sort(header.events.POS);
    header.events.POS = header.events.POS(tmp);
    header.events.DUR = header.events.DUR(tmp);
    header.events.TYP = header.events.TYP(tmp);
    header.events.OFF = header.events.OFF(tmp);
    clearvars tmp bad settingsB duration
    domarkers = true;
else
    domarkers = false;
end

if isfield(cfg,'Output_file') & exist(fullfile(cfg.Output_filepath,cfg.Output_file),'file')
    [~,~,~,Output_fileS] = lab_filename(cfg.Output_file);
    Output_fileI = fullfile(cfg.Output_filepath,Output_fileS);
    Output_file = fullfile(cfg.Output_filepath,cfg.Output_file);
elseif isfield(cfg,'EEG_file') & exist(fullfile(cfg.EEG_filepath,cfg.EEG_file),'file')
    [~,~,~,EEG_fileS] = lab_filename(cfg.EEG_file);
    Output_fileI = fullfile(cfg.EEG_filepath,EEG_fileS);
    Output_file = fullfile(cfg.EEG_filepath,cfg.EEG_file);
else
    Output_fileI = '';
end
if ~isempty(Output_fileI) & exist(Output_fileI,'file')
    headertmp = lab_read_eeginfo(Output_fileI);
    if isfield(headertmp,'numdatachannels') & headertmp.numdatachannels == header.numdatachannels
        headertmp.badchans = header.badchans;
        lab_write_eeginfo(Output_fileI,headertmp);
        if domarkers == true
            lab_write_mrk([Output_file '.mrk'],header);
        end
    end
end

return