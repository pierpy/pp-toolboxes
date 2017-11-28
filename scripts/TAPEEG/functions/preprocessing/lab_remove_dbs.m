function [data,header,cfg,skipprocessing] = lab_remove_dbs(data,header,cfg)
    
disp('Remove DBS-stimulator artifacts')
    
skipprocessing = 0;

if ~exist('cfg','var')
    cfg = [];
end
if ~exist('header','var')
    header = lab_create_header(data);
end

if~isfield(cfg,'STIM') | isempty(cfg.STIM)
    [cfg,skipprocessing] = lab_set_remove_dbs(cfg,header);
    if skipprocessing == 1
        return
    end
end

if ~isfield(header,'numdatachannels')
    header.numdatachannels = header.numchannels;
end

if ~isfield(cfg,'EEG_filepath')
    if isfield(header,'EEG_filepath')
        cfg.EEG_filepath = header.EEG_filepath;
    else
        cfg.EEG_filepath = pwd;
    end
end
if ~isfield(cfg,'EEG_file')
    if isfield(header,'EEG_file')
        cfg.EEG_file = header.EEG_file;
    else
        cfg.EEG_file = 'File.sef';
    end
end

% highpass filter
disp('   high pass filter')
if isfield(cfg,'FILT') & isfield(cfg.FILT,'Filter')
    settings.FILT = cfg.FILT;
    for i = 1:length(settings.FILT.Filter)
        settings.FILT.Filter(i).lowpass = [];
        settings.FILT.Filter(i).notch = [];
    end
else
    settings.FILT.filtermode = 'butter';
    settings.FILT.filtorder = 2;
    settings.FILT.Filter.firstchan = 1;
    settings.FILT.Filter.lastchan = header.numdatachannels;
    settings.FILT.Filter.highpass = 1;
    settings.FILT.Filter.lowpass = [];
    settings.FILT.Filter.notch = [];
end
[data,header] = lab_filtering(data,header,settings);

% find stimulator artifacts
disp('   find dbs artifacts')
if isfield(cfg.STIM,'DETECT') & ~isempty(cfg.STIM.DETECT)
    [~,headerS,cfg,skipprocessing] = lab_detect_dbsstim(data,header,cfg);
    if skipprocessing == 1
        disp('   no dbs artifact found, skip ICA')
        return
    end
end

% do principal component analysis
if isfield(cfg.STIM,'PCA') & cfg.STIM.PCA == true
    if mean(mean(abs(data))) > 400
        disp('   remove principal component')
        Idx = setdiff(1:header.numdatachannels,find(max(abs(data),[],2) == 0));
        [weights, sphere, ~, ~, ~, ~, activations] = runica(data(Idx,:),'verbose','on','extended',0,'pca',1);
        W = pinv(weights*sphere);
        data(Idx,:) = data(Idx,:) - W*activations;
    else
        disp('   skip PCA, not enough power in data')
    end
end

if isfield(cfg.STIM,'ICA') & ~isempty(cfg.STIM.ICA)
    % do ICA and detect components by averaging on stimulator artifacts
    disp('   do ICA')
    if isfield(cfg,'exclude')
        cfg.STIM.exclude = cfg.exclude;
    end
    [~,~,format,cfg.EEG_fileS] = lab_filename(cfg.EEG_file);
    cfg.STIM.EEG_file = [cfg.EEG_fileS '_S' num2str(i) '.' format];
    cfg.STIM.EEG_filepath = fullfile(cfg.EEG_filepath,cfg.STIM.folder);
    [~,headerICA,cfg.STIM] = lab_ICAstart(data,headerS,cfg.STIM);
    
    % Remove components with stimulator and re-calculate EEG
    Filename = fullfile(cfg.STIM.ICA_filepath,[cfg.STIM.ICA_file(1:end-3) 'mat']);
    cfg.STIM.ICABACK.FILEBAD = true;
    cfg.STIM.ICABACK.BAD = [];
    cfg.STIM.ICABACK.ACTIVATIONS = [];
    cfg.STIM.ICABACK.dobacktransform = true;
    cfg.STIM.ICABACK.foldername = 'ICAresult';
    cfg.STIM.ICABACK.IS = [];
    cfg.STIM.ICABACK.PREPROCESSING  = [];
    [data,~,cfg.STIM] = lab_ICAback(Filename,cfg.STIM,headerICA.badchans);
end

% lowpass filter
disp('   low pass filter')
if isfield(cfg,'FILT') & isfield(cfg.FILT,'Filter')
    settings.FILT = cfg.FILT;
    for i = 1:length(settings.FILT.Filter)
        settings.FILT.Filter(i).highpass = [];
    end
else
    settings.FILT.filtermode = 'firls';
    settings.FILT.filtorder = 0;
    settings.FILT.Filter.firstchan = 1;
    settings.FILT.Filter.lastchan = header.numdatachannels;
    settings.FILT.Filter.highpass = [];
    settings.FILT.Filter.lowpass = 70;
    settings.FILT.Filter.notch = [];
end
[data,header] = lab_filtering(data,header,settings);

end