% Match LOCS and channels in EEG/MEG-file
%
% [data,header] = lab_match_locs_eeg(data,header,cfg)
%
% written by F. Hatz 2014

function [data,header] = lab_match_locs_eeg(data,header,cfg)

if ~exist('cfg','var') | ~isfield(cfg,'LOCS')
    cfg.LOCS = [];
end

% Get filename and path
if isfield(cfg,'EEG_file')
    filename = cfg.EEG_file;
    filepath = cfg.EEG_filepath;
elseif isfield(header,'EEG_file')
    filename = header.EEG_file;
    filepath = header.EEG_filepath;
else
    filename = 'EEG';
    filepath = pwd;
end
[~,~,~,filenameS] = lab_filename(filename);

if isfield(cfg.LOCS,'locs') & ~isempty(cfg.LOCS.locs)
    LOCS = cfg.LOCS.locs;
elseif ~isfield(cfg.LOCS,'filelocs') | cfg.LOCS.filelocs == true
    if exist(fullfile(filepath,[filenameS '.els']),'file')
        LOCS = lab_read_locs(fullfile(filepath,[filenameS '.els']));
    elseif exist(fullfile(filepath,[filenameS '.elc']),'file')
        LOCS = lab_read_locs(fullfile(filepath,[filenameS '.elc']));
    elseif exist(fullfile(filepath,[filenameS '.xyz']),'file')
        LOCS = lab_read_locs(fullfile(filepath,[filenameS '.xyz']));
    elseif exist(fullfile(filepath,[filenameS '.sfp']),'file')
        LOCS = lab_read_locs(fullfile(filepath,[filenameS '.sfp']));
    elseif exist(fullfile(filepath,[filenameS '.spi']),'file')
        LOCS = lab_read_locs(fullfile(filepath,[filenameS '.spi']));
    elseif exist(fullfile(filepath,'electrodes.els'),'file')
        LOCS = lab_read_locs(fullfile(filepath,'electrodes.els'));
    elseif exist(fullfile(filepath,'electrodes.elc'),'file')
        LOCS = lab_read_locs(fullfile(filepath,'electrodes.elc'));
    elseif exist(fullfile(filepath,'electrodes.sfp'),'file')
        LOCS = lab_read_locs(fullfile(filepath,'electrodes.sfp'));
    elseif exist(fullfile(filepath,'electrodes.xyz'),'file')
        LOCS = lab_read_locs(fullfile(filepath,'electrodes.xyz'));
    else
        LOCS = [];
    end
elseif isfield(header,'locs') & ~isempty(header.locs)
    LOCS = header.locs;
else
    LOCS = [];
end
if isempty(LOCS)
    LOCS = lab_read_locs;
end

if ~isfield(LOCS,'labels') | ~isfield(header,'channels')
    return
end

disp('   Match channels and locs')

channels = textscan(header.channels','%s');
channels = channels{1,1};
labels = LOCS.labels;
channels = lab_correct_labels(channels,labels);

Ilabels = [];
Ichannels = [];
chanlabels = {};
for i = 1:length(labels)
    tmp = find(~cellfun(@isempty,strfind(upper(channels),upper(labels{i}))));
    tmp = setdiff(tmp,Ichannels);
    if ~isempty(tmp)
        if length(tmp) > 1
            tmp2 = find(strcmpi(channels(tmp),labels{i}));
            if ~isempty(tmp2)
                tmp = tmp2(1);
            else
                tmp = tmp(1);
            end
        end
        Ilabels = [Ilabels i]; %#ok<AGROW>
        chanlabels{end+1,1} = labels{i}; %#ok<AGROW>
        Ichannels = [Ichannels tmp]; %#ok<AGROW>
    end
end
if isempty(Ichannels)
    data = [];
    header = [];
    disp('    Abort: loc-file and input-file not matching');
    return
end

Iaux = setdiff(1:length(channels),Ichannels);
auxlabels = channels(Iaux);
auxlabels = auxlabels(:);
Iaux = Iaux(:)';
if ~isempty(auxlabels)
    tmp = find(~cellfun(@isempty,strfind(upper(auxlabels),'EOG')));
    if ~isempty(tmp)
        if length(tmp) > 1
            tmp2 = find(~cellfun(@isempty,strfind(upper(auxlabels(tmp)),'VEOG')));
            if ~isempty(tmp2)
                tmp = tmp(tmp2(1));
            end
        end
        if length(auxlabels) > 1
            Iaux = [Iaux(tmp(1)) setdiff(Iaux,Iaux(tmp(1)),'stable')];
            auxlabels = cat(1,cellstr('EOG'),setdiff(auxlabels,auxlabels(tmp(1)),'stable'));
        end
        header.eog_ch = length(Ichannels) + 1;
    end
    tmp = find(~cellfun(@isempty,strfind(upper(auxlabels),'ECG')));
    if ~isempty(tmp)
        if length(auxlabels) > 1
            Iaux = [Iaux(tmp) setdiff(Iaux,Iaux(tmp),'stable')];
            auxlabels = cat(1,cellstr('ECG'),setdiff(auxlabels,auxlabels(tmp),'stable'));
        end
        header.ecg_ch = length(Ichannels) + 1;
        if isfield(header,'eog_ch')
            header.eog_ch = header.eog_ch + 1;
        end
    end
    tmp = find(~cellfun(@isempty,strfind(upper(auxlabels),'EKG')));
    if ~isempty(tmp)
        if length(auxlabels) > 1
            Iaux = [Iaux(tmp) setdiff(Iaux,Iaux(tmp),'stable')];
            auxlabels = cat(1,cellstr('ECG'),setdiff(auxlabels,auxlabels(tmp),'stable'));
        end
        header.ecg_ch = length(Ichannels) + 1;
        if isfield(header,'eog_ch')
            header.eog_ch = header.eog_ch + 1;
        end
    end
end
MATCHLOCS = LOCS;
LOCS.labels = LOCS.labels(1,Ilabels);
LOCS.x = LOCS.x(1,Ilabels);
LOCS.y = LOCS.y(1,Ilabels);
LOCS.z = LOCS.z(1,Ilabels);
LOCS = lab_locs2sph(LOCS);
if ~isempty(Iaux)
    LOCS.aux = length(Iaux);
    LOCS.auxlabels = auxlabels;
end
LOCS.MATCHLOCS = MATCHLOCS;
header.locs = LOCS;
header.numauxchannels = length(Iaux);
header.numdatachannels = length(Ichannels);
header.channels = char(cat(1,chanlabels,auxlabels));

if ~isempty(data)
    data = data([Ichannels Iaux],:);
end
header.numchannels = size(data,1);

end