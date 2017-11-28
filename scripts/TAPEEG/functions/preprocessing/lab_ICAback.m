% Calculate resulting EEG after ICA analysis
%
% [data,header,cfg] = lab_ICAback(filename,activations_exclude,cfgorig)
%
% filename            = .mat file stored with preprocessing ICA
% activations         = matrix (activations x timeframes)
% activations_exclude = index of activations to exclude / auto = automatic search
% cfgorig             = structure with config (optional)
% dosearch            = if 1 -> only bad activations search is done and
%                               written to file 'exclude~.txt'
%                       if 2 -> bad activations search and calculation of
%                               resulting EEG
%
% written by F. Hatz 2012


function [data,header,cfg,skipprocessing] = lab_ICAback(filename,cfgorig,activations_exclude)

skipprocessing = 0;

disp('ICA backtransform')

if ~exist('filename','var') | isempty(filename)
    [filename,filepath] = uigetfile('*.mat','Select ICA.mat-file');
    filename = fullfile(filepath,filename);
    clearvars filepath
end

% Load ICA-result-file
load(filename);
if ~exist('header','var') | ~exist('W','var')
    if exist('cfgorig','var') & isfield(cfgorig,'EEG_file')
        disp(['   No valid ICA-result-file (' cfgorig.EEG_file ')'])
    else
        disp('   No valid ICA-result-file')
    end
    data = [];
    header = [];
    skipprocessing = 1;
    return
end

if ~exist('cfgorig','var') | ~isfield(cfgorig,'ICABACK') | isempty(cfgorig.ICABACK)
    [cfgorig,skipprocessing] = lab_set_ICAback([],header); %#ok<NODEF>
    pause(0.2);
    if skipprocessing == 1 | ~isfield(cfgorig,'ICABACK')
        data = [];
        header = [];
        cfg = [];
        skipprocessing = 1;
        return
    end
end
ICABACK = cfgorig.ICABACK;
if isempty(ICABACK) | ~isfield(ICABACK,'dobacktransform')
    data = [];
    header = [];
    cfg = [];
    skipprocessing = 1;
    return
end

if exist('cfgorig','var') & isfield(cfgorig,'ICA') & isfield(cfgorig.ICA,'good')
    cfg.ICA.good = cfgorig.ICA.good;
else
    cfg.ICA.good = 1;
end

% Correct for old data
if exist('ICA_file','var')
    cfg.ICA_file = ICA_file; %#ok<NODEF>
    cfg.exclude = exclude;
    clearvars ICA_file exclude
end
if exist('good','var')
    header.goodchans = good;
    header.badchans = bad;
    header.bad.All = bad;
    clearvars good bad
end
if ~isfield (header,'ecg_ch')
    header.ecg_ch = 0;
end
if ~isfield(cfg,'EEG_file')
    cfg.EEG_file = [cfg.ICA_file(1:end-8) '.sef'];
end
if isfield(header,'badchans') & max(header.badchans == 0);
    if ~isempty(find(header.badchans ~= 0,1))
        header.badchans = header.badchans(header.badchans ~= 0);
    else
        header.badchans = [];
    end
end

% Correct file_paths
[~,cfg.ICA_filepath] = lab_filename(filename);
tmp=strfind(cfg.ICA_filepath,filesep);
cfg.settings_path = cfg.ICA_filepath(1:tmp(end-2));
cfg.EEG_filepath = cfg.ICA_filepath(1:tmp(end-1));
clearvars tmp pathtmp;

if ~exist(fullfile(cfg.ICA_filepath,cfg.ICA_file),'file') & exist([filename(1:end-3) 'sef'],'file')
    [~,~,~,cfg.ICA_fileS] = lab_filename(filename);
    cfg.ICA_file = [cfg.ICA_fileS '.sef'];
    clearvars tmp
end
if isfield(cfg,'EEG_output_file')
    if ~exist(fullfile(cfg.EEG_filepath,cfg.EEG_output_file),'file') & exist(fullfile(cfg.EEG_filepath,[cfg.ICA_file(1:end-8) 'filt.sef']),'file')
        cfg.EEG_output_file = [cfg.ICA_file(1:end-8) 'filt.sef'];
    end
else
    if ~exist(fullfile(cfg.EEG_filepath,cfg.EEG_file),'file') & exist(fullfile(cfg.EEG_filepath,[cfg.ICA_file(1:end-8) 'filt.sef']),'file')
        cfg.EEG_file = [cfg.ICA_file(1:end-8) 'filt.sef'];
    end
end

if isfield(cfg,'EEG_file')
    header.EEG_file = cfg.EEG_file;
    header.EEG_filepath = cfg.EEG_filepath;
end

% Read electrodes-file when missing
if ~isfield(header,'locs') & exist(fullfile(cfg.settings_path,'electrodes.sfp'),'file')
    header.locs = lab_read_locs(fullfile(cfg.settings_path,'electrodes.sfp'));
end

%--------------------------------------------------------------------------
% LOAD ICA COMPONENTS
%--------------------------------------------------------------------------
disp('   Load activations');
[activations,headeract] = lab_read_data(fullfile(cfg.ICA_filepath,cfg.ICA_file));
if isempty(activations)
    disp('   No valid ICA-activations-file')
    data = [];
    header = [];
    cfg = cfgorig;
    skipprocessing = 1;
    return
end
header = lab_mix_markers(header,headeract);

%--------------------------------------------------------------------------
% Set Output folder
%--------------------------------------------------------------------------
if isfield(ICABACK,'foldername') & ~isempty(ICABACK.foldername)
    cfg.finalfolder = ICABACK.foldername;
else
    cfg.finalfolder = 'ICAresult';
end

%--------------------------------------------------------------------------
% Read original EEG for Bad/Aux-channels
%--------------------------------------------------------------------------
if header.numchannels > size(W,1) %#ok<NODEF>
    disp('   Read original EEG for Bad/Aux-channels');
    if ~isfield(cfg,'EEG_output_file')
        cfg.EEG_output_file = [cfg.EEG_file(1:end-4) 'filt.sef'];
    end
    if exist(fullfile(cfg.ICA_filepath,cfg.EEG_output_file),'file')
        EEG_filt_file = fullfile(cfg.ICA_filepath,cfg.EEG_output_file);
    elseif exist(fullfile(cfg.EEG_filepath,cfg.EEG_output_file),'file')
        EEG_filt_file = fullfile(cfg.EEG_filepath,cfg.EEG_output_file);
    end
    [data,headerfilt] = lab_read_data(EEG_filt_file);
    if isempty(data) | size(data,1) ~= header.numchannels
        disp('   Original EEG for Bad/Aux-channels not found, replace by zeros');
        data = zeros(header.numchannels,header.numtimeframes);
        EEG_filt_file = '';
    end
else
    data = zeros(header.numchannels,header.numtimeframes);
    EEG_filt_file = '';
end

%--------------------------------------------------------------------------
% Define activations to exclude
%--------------------------------------------------------------------------
disp('   Define activations to exclude')
if ~exist('activations_exclude','var')
    if isfield(ICABACK,'ACTIVATIONS') & isempty(ICABACK.ACTIVATIONS) & ...
            isfield(ICABACK,'FILEBAD') & ICABACK.FILEBAD == false & ...
            isfield(ICABACK,'BAD') & isempty(ICABACK.BAD)
        disp('     no activations excluded')
        activations_exclude = [];
        flagwrite = false;
    elseif isfield(ICABACK,'ACTIVATIONS') & ~isempty(ICABACK.ACTIVATIONS)
        disp('     fixed bad activations (configuration)')
        activations_exclude = ICABACK.ACTIVATIONS;
        headeract.badchans = activations_exclude;
        headeract.goodchans = setdiff(1:size(activations,1),headeract.badchans);
        headeract.goodchans = headeract.goodchans(:)';
        flagwrite = true;
    elseif isfield(ICABACK,'BAD') & ~isempty(ICABACK.BAD)
        disp('     detect bad activations')
        if isfield(ICABACK.BAD,'ecg_ch') & ICABACK.BAD.ecg_ch > 0 & ICABACK.BAD.ecg_ch <= size(data,1)
            ICABACK.BAD.ecg_chan = data(ICABACK.BAD.ecg_ch,:);
        end
        if isfield(ICABACK.BAD,'eog') & ~isempty(ICABACK.BAD.eog)
            ICABACK.BAD = lab_calculate_eog(data,header,ICABACK.BAD);
        end
        headeract.W = W;
        if isfield(headeract,'locs')
            headeract = rmfield(headeract,'locs');
        end
        EEG_file = cfg.EEG_file;
        EEG_filepath = cfg.EEG_filepath;
        cfg.EEG_file = cfg.ICA_file;
        cfg.EEG_filepath = cfg.ICA_filepath;
        [activations_exclude,~,ICABACK.BAD] = lab_detect_bad(activations,headeract,ICABACK.BAD,cfg);
        headeract.badchans = activations_exclude;
        headeract.goodchans = setdiff(1:size(activations,1),headeract.badchans);
        headeract.goodchans = headeract.goodchans(:)';
        cfg.EEG_file = EEG_file;
        cfg.EEG_filepath = EEG_filepath;
        flagwrite = true;
    elseif isfield(headeract,'badchans') & (~isfield(ICABACK,'FILEBAD') | ICABACK.FILEBAD == true)
        disp('     bad activations by file')
        activations_exclude = headeract.badchans;
        flagwrite = false;
    elseif ICABACK.dobacktransform == true & isfield(cfgorig,'MAIN') & ...
            isfield(cfg.MAIN,'auto') & cfg.MAIN.auto == 1
        disp('   Abort: no valid data for excluded activations')
        data = [];
        header = [];
        cfg = cfgorig;
        skipprocessing = 1;
        return
    elseif ICABACK.dobacktransform == true
        disp ('   Ask for activations to exclude')
        defaultanswer = {num2str(0)};
        answer = inputdlg({'Excluded activations'},'Excluded activations',[1 50],defaultanswer);
        pause(0.1)
        if size(answer,1) > 0
            activations_exclude = str2num(answer{1,1}); %#ok<ST2NM>
        else
            disp('   Abort: wrong data for excluded activations')
            data = [];
            header = [];
            cfg = cfgorig;
            skipprocessing = 1;
            return
        end
        clearvars defaultanswer answer
        flagwrite = true;
    else
        disp('   Abort: wrong data for processing')
        data = [];
        header = [];
        cfg = cfgorig;
        skipprocessing = 1;
        return
    end
    if flagwrite == true
        % write bad activations to files
        if exist(fullfile(cfg.ICA_filepath,[cfg.ICA_file(1:end-4) '_exclude.txt']),'file')
            delete(fullfile(cfg.ICA_filepath,[cfg.ICA_file(1:end-4) '_exclude.txt']));
        end
        disp('    write exclude~.txt')
        fid=fopen(fullfile(cfg.ICA_filepath,[cfg.ICA_file(1:end-4) '_exclude~.txt']),'w');
        fprintf(fid,num2str(activations_exclude));
        fclose(fid);
        disp('    write .info')
        lab_write_eeginfo(fullfile(cfg.ICA_filepath,[cfg.ICA_file(1:end-4) '.sef']),headeract);
    end
    clearvars flagwrite
end
if ~isnumeric(activations_exclude)
    disp('   Abort: wrong data for excluded activations')
    data = [];
    header = [];
    cfg = cfgorig;
    skipprocessing = 1;
    return
end
activations_exclude = activations_exclude(activations_exclude>0);
if max(activations_exclude) > size(W,1)
    disp(['   Abort: max activation to exclude (' num2str(max(activations_exclude)) ') to large'])
    data = [];
    header = [];
    cfg = cfgorig;
    skipprocessing = 1;
    return
end

%--------------------------------------------------------------------------
% Skip processing if backtransform disabled
%--------------------------------------------------------------------------
if ICABACK.dobacktransform == false
    data = [];
    header = [];
    cfg = cfgorig;
    skipprocessing = 1;
    return
end

%--------------------------------------------------------------------------
% Remove selected ICA components
%--------------------------------------------------------------------------
if ~exist('ICAchans','var')
    if cfg.ICA.good == 1 & isfield(header,'goodchans');
        ICAchans = header.goodchans;
    else
        ICAchans = 1:header.numdatachannels;
        if isfield(header,'ref_chan') & isnumeric(header.ref_chan)
            ICAchans = setdiff(ICAchans,header.ref_chan);
        end
        if isfield(cfg,'exclude') & ~isempty(cfg.exclude)
            ICAchans = setdiff(ICAchans,cfg.exclude);
        end
    end
end
if length(ICAchans) ~= size(W,1)
    if size(W,1) == header.numdatachannels
        ICAchans = 1:header.numdatachannels;
    else
        disp('   Abort, mismatch in data')
        data = [];
        header = [];
        cfg = cfgorig;
        skipprocessing = 1;
        return
    end
end
disp('   Calculate EEG');
if min(activations_exclude) > 0
    disp(['     activations excluded: ' num2str(activations_exclude(:)')])
    W(:,activations_exclude)=0;
else
    disp('     activations excluded: none')
end
if isreal(W)
    data(ICAchans,:)=W*activations;
    header.activationsexcluded = length(activations_exclude);
else
    disp('   Abort, ICA failed by giving complex results')
    data = [];
    header = [];
    cfg = cfgorig;
    skipprocessing = 1;
    return
end

%--------------------------------------------------------------------------
% Take bad channels info from *filt
%--------------------------------------------------------------------------
if exist('headerfilt','var') & isfield(headerfilt,'badchans')
    header.badchans = headerfilt.badchans;
    header.goodchans = setdiff(header.goodchans,headerfilt.badchans);
    header.goodchans = header.goodchans(:)';
end

%--------------------------------------------------------------------------
% Define Output_folder
%--------------------------------------------------------------------------
cfg.Output_filepath = fullfile(cfg.ICA_filepath,cfg.finalfolder);
warning off %#ok<WNOFF>
mkdir (cfg.Output_filepath);
warning on %#ok<WNON>
cfg.Output_file=[cfg.ICA_file(1:end-8) '.sef'];

%--------------------------------------------------------------------------
% Finalize cfg
%--------------------------------------------------------------------------
ICA_file = cfg.ICA_file;
cfgorig.EEG_file = cfg.EEG_file;
cfgorig.EEG_filepath = cfg.EEG_filepath;
cfgorig.Output_file = cfg.Output_file;
cfgorig.Output_filepath = cfg.Output_filepath;
cfg = cfgorig;

%--------------------------------------------------------------------------
% Do dipol fit, if selected
%--------------------------------------------------------------------------
if isfield(ICABACK,'IS') & ~isempty(ICABACK.IS)
    include = setdiff(1:size(activations,1),activations_exclude);
    W = W(:,include);
    headeract.W = W;
    [activations,headeract] = lab_reduce_channels(activations,headeract,include,true);
    lab_ICA_dipolfit(activations,headeract,cfg);
end
clearvars W;

%--------------------------------------------------------------------------
% Verbose file
%--------------------------------------------------------------------------
[~,~,~,Verbose_file] = lab_filename(cfg.Output_file);
Verbose_file = fullfile(cfg.Output_filepath,[Verbose_file '_ICAback.vrb']);
fid=fopen(Verbose_file,'w');
fprintf(fid,'ICA backtransform\n');
fprintf(fid,datestr(now,0));
fprintf(fid,'\n');
fprintf(fid,['EEG input file: ' lab_filename(EEG_filt_file)]);
fprintf(fid,'\n\n');
fprintf(fid,['Activations input file: ' ICA_file]);
fprintf(fid,'\n\n');
if isempty(activations_exclude)
    fprintf(fid,'Excluded activations: none');
else
    fprintf(fid,['Excluded activations: ' num2str(activations_exclude) '\n']);
    activations_included = setdiff(1:size(activations,1),activations_exclude);
    activations_included = activations_included(:)';
    fprintf(fid,['Included activations: ' num2str(activations_included) '\n']);
end
fclose(fid);
clearvars activations

% do preprocessing if enabled
if isfield(cfg.ICABACK,'PREPROCESSING') & ~isempty(cfg.ICABACK.PREPROCESSING)
    try
        cfg.ICABACK.PREPROCESSING.EEG_file = cfg.EEG_file;
        cfg.ICABACK.PREPROCESSING.EEG_filepath = cfg.EEG_filepath;
        cfg.ICABACK.PREPROCESSING.Output_file = cfg.Output_file;
        cfg.ICABACK.PREPROCESSING.Output_filepath = cfg.Output_filepath;
        cfg.ICABACK.PREPROCESSING.settings_path = cfg.settings_path;
        if isfield(cfg,'MAIN') & isfield(cfg.MAIN,'auto')
            cfg.ICABACK.PREPROCESSING.MAIN.auto = cfg.MAIN.auto;
        else
            cfg.ICABACK.PREPROCESSING.MAIN.auto = 0;
        end
        [data,header,cfg.ICABACK.PREPROCESSING,skipprocessing] = lab_preprocessing(data,header,cfg.ICABACK.PREPROCESSING);
        cfg.ICABACK.PREPROCESSING = rmfield(cfg.ICABACK.PREPROCESSING,'EEG_file');
        cfg.ICABACK.PREPROCESSING = rmfield(cfg.ICABACK.PREPROCESSING,'EEG_filepath');
        cfg.ICABACK.PREPROCESSING = rmfield(cfg.ICABACK.PREPROCESSING,'Output_file');
        cfg.ICABACK.PREPROCESSING = rmfield(cfg.ICABACK.PREPROCESSING,'Output_filepath');
        cfg.ICABACK.PREPROCESSING = rmfield(cfg.ICABACK.PREPROCESSING,'settings_path');
        cfg.ICABACK.PREPROCESSING = rmfield(cfg.ICABACK.PREPROCESSING,'MAIN');
    catch err
        disp(getReport(err))
        skipprocessing = 1;
    end
end

return

