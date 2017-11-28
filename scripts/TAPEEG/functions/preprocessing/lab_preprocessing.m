% Preprocessing of eeg/meg-data, script called by lab_tapeeg
%
% [data,header,cfg,skipprocessing,skiperror] = lab_preprocessing(data,header,cfg,calc,skipfilt)
%
% Written by F. Hatz 2013 Neurology Basel

function [data,header,cfg,skipprocessing,skiperror] = lab_preprocessing(data,header,cfg,calc,skipfilt)
disp ('-Preprocessing-')

skipprocessing = 0;
skiperror = 0;
dodbs = false;

if ~exist('skipfilt','var')
    skipfilt = 0;
end
if ~exist('calc','var')
    calc = [];
end
if ~exist('cfg','var')
    cfg = [];
end

if ~isfield(cfg,'EEG_file')
    cfg.EEG_file = header.EEG_file;
    cfg.EEG_filepath = header.EEG_filepath;
end

if skipfilt == 0
    % Scaling
    if isfield(cfg,'SCALE') & ~isempty(cfg.SCALE)
        [data,header,cfg,skipprocessing] = lab_scale_data(data,header,cfg);
        if skipprocessing == 1
            skiperror = 1;
            return
        end
    end
    
    % Remove DBS-stimulator
    if isfield(cfg,'STIM') & ~isempty(cfg.STIM)
        dodbs = true;
        [data,header,cfg,skipprocessing] = lab_remove_dbs(data,header,cfg);
        if skipprocessing == 1
            skiperror = 1;
            return
        end
    end
    
    % Filtering
    if isfield(cfg,'FILT') & ~isempty(cfg.FILT) & dodbs == false
        [data,header,cfg,skipprocessing] = lab_filtering(data,header,cfg);
        if skipprocessing == 1
            skiperror = 1;
            return
        end
    end
    
    % Resample
    if isfield(cfg,'RESAMPLE') & isfield(cfg.RESAMPLE,'new_samplingrate') & ~isempty(cfg.RESAMPLE.new_samplingrate)
        [~,~,~,verbosefile] = lab_filename(cfg.EEG_file);
        verbosefile = fullfile(cfg.EEG_filepath,[verbosefile '_resample.vrb']);
        [data,header,cfg,skipprocessing] = lab_resample_data(data,header,cfg,verbosefile);
        if skipprocessing == 1
            skiperror = 1;
            return
        end
    end
    
    % Stitching
    if isfield(cfg,'STITCH') & ~isempty(cfg.STITCH)
        [data,header,cfg,skipprocessing] = lab_stitching(data,header,cfg,calc);
        if skipprocessing == 1
            return
        end
    end
    
    % find good/bad channels
    if isfield(cfg,'BADELEC') & ~isempty(cfg.BADELEC)
        [header,cfg] = lab_find_good(data,header,cfg);
        if ~isfield(header,'bad') | ~isfield(header.bad,'error') | header.bad.error == 1
            skipprocessing = 1;
            skiperror = 1;
            return
        end
    elseif ~isfield(header,'goodchans')
        if isfield(header,'numdatachannels')
            header.goodchans = 1:header.numdatachannels;
        else
            header.goodchans = 1:header.numchannels;
        end
        header.badchans = [];
    end
end

% ICA
if isfield(cfg,'ICA') & ~isempty(cfg.ICA)
    [~,headerICA,cfg] = lab_ICAstart(data,header,cfg);
    if isempty(headerICA)
        skiperror = 1;
    end
    if isfield(cfg.ICA,'ICABACK') & ~isempty(cfg.ICA.ICABACK)
        if isempty(headerICA)
            disp('ICA failed, skip data')
            data = [];
            header = [];
            skipprocessing = 1;
        elseif isfield(headerICA,'W') & isfield(headerICA,'badchans') & isfield(headerICA,'bad') & ...
                isfield(headerICA.bad,'error') & headerICA.bad.error == 0
            Filename = fullfile(cfg.ICA_filepath,[cfg.ICA_file(1:end-3) 'mat']);
            cfg2 = cfg;
            cfg2.ICABACK = cfg.ICA.ICABACK;
            [data,header,cfg2,skipprocessing] = lab_ICAback(Filename,cfg2,headerICA.badchans);
            cfg.ICA.ICABACK = cfg2.ICABACK;
            clearvars cfg2 Filename
        else
            disp('ICA: Detection of bad activations failed, please select manuallly')
            clearvars W activations badAll;
            skipprocessing = 1;
            skiperror = 1;
        end
    elseif ~isempty(data)
        disp('ICA: Manually select activations to exclude and restart...')
        cfg.ICA.automated = false;
        clearvars W activations badAll;
        skipprocessing = 1;
        return
    end
end

% Add noise
if isfield(cfg,'NOISE') & ~isempty(cfg.NOISE) & skipprocessing == 0 & ~isempty(data)
    [data,header,cfg,skipprocessing] = lab_add_noise(data,header,cfg);
    if skipprocessing == 1
        skiperror = 1;
        return
    end
end

% Smoothing
if isfield(cfg,'SMOOTH') & ~isempty(cfg.SMOOTH) & skipprocessing == 0 & ~isempty(data)
    [data,header,cfg,skipprocessing] = lab_smoothing(data,header,cfg,true);
    if skipprocessing == 1
        skiperror = 1;
        return
    end
end

% Stitch All
if isfield(cfg,'STITCHALL') & ~isempty(cfg.STITCHALL)
    [data,header,cfg,skipprocessing] = lab_stitchingall(data,header,cfg);
end
if skipprocessing == 1
    return
end

% Interpolate bad
if isfield(cfg,'INTERPOLATE') & ~isempty(cfg.INTERPOLATE) & ...
        isfield(cfg.INTERPOLATE,'method') & ~isempty(data)
     disp ('   Interpolate bad channels')
    [data,header,cfg.INTERPOLATE.method] = lab_interpolate_bad(data,header,cfg.INTERPOLATE.method);
end

% Detect eyeblinks
if isfield(cfg,'EYEBLINKS') & ~isempty(cfg.EYEBLINKS) & ...
        isfield(cfg.EYEBLINKS,'eog') & ~isempty(data)
    disp ('   Detect Eyeblinks')
    [data,header,cfg] = lab_detect_eyeblinks(data,header,cfg);
end

% Individual Frequency Bands
if isfield(cfg,'IFREQ') & ~isempty(cfg.IFREQ) & ~isempty(data)
    disp ('   Detect Individual Frequency Bands')
    [data,header,cfg] = lab_indiv_freqbands(data,header,cfg);
end


return