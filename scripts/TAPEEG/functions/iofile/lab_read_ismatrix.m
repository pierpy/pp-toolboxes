function settings = lab_read_ismatrix(settings,cfg,numdatachannels)

if isempty(settings.IS_file) | ~exist(settings.IS_file,'file')
    return
end

% read IS matrix
if strcmp(settings.IS_file(end-2:end),'.is')
    if ~isfield(settings,'matrix') & exist(settings.IS_file,'file')
        disp('   Read Cartool-IS-matrix')
        if ~isfield(settings,'regularization') | isempty(settings.regularization)
            settings.regularization = 6;
        end
        [settings.matrix,settings.regularization] = lab_read_is(settings.IS_file,settings.regularization);
    end
    settings.LF = [];
elseif strcmp(settings.IS_file(end-5:end),'.spinv')
    if ~isfield(settings,'matrix') & exist(settings.IS_file,'file')
        disp('   Read sLoreta-matrix')
        [settings.matrix,settings.regularization] = lab_read_is(settings.IS_file);
    end
    settings.LF = [];
elseif strcmp(settings.IS_file(end-3:end),'.bin') & isfield(settings,'slocs') & exist(settings.IS_file,'file')
    disp('   Read Leadfield (.bin)')
    if exist('numdatachannels','var') & ~isempty(numdatachannels) & isnumeric(numdatachannels)
        numchannels = numdatachannels;
    elseif isfield(cfg,'numdatachans') & ~isempty(cfg.numdatachans)
        numchannels = cfg.numdatachans;
    elseif isfield(cfg,'EXTRA') & isfield(cfg.EXTRA,'numdatachans') & ~isempty(cfg.EXTRA.numdatachans)
        numchannels = cfg.EXTRA.numdatachans;
    else
        disp ('    Ask for number of data channels')
        answer=inputdlg({'Number of data channels'},'data channels',[1 40]);
        pause(0.1)
        if size(answer) > 0
            numchannels = str2num(answer{1,1});
        else
            return
        end
        clearvars answer name numlines prompt
    end
    if isstruct(settings.slocs)
        numsp = size(settings.slocs.x,2);
    else
        numsp = size(settings.slocs,2);
    end
    try 
        settings.LF.LFbin = lab_read_LFbin(settings.IS_file,numchannels,numsp);
    catch %#ok<CTCH>
        disp('    Abort: error reading Leadfield')
        settings.LF.LFbin = [];
    end
elseif strcmp(settings.IS_file(end-3:end),'.mat') & isfield(settings,'slocs') & exist(settings.IS_file,'file')
    disp('   Read Headmodel (.mat)')
    load(settings.IS_file);
    if exist('vol','var')
        settings.LF.LFbin = [];
        settings.LF.vol = vol;
    end
elseif isfield(settings,'LF') & isfield(settings.LF,'LFbin')
    settings.LF.LFbin = [];
end