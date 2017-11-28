% TAPEEG start-script
%
% Parameters: 'tapeeg'      = start main application to process
%                             files or folders
%             'statistics'  = statistics, input are xls-files
%             'plot'        = plot xls-files or results of permutation
%                             in singal or source space
%   'folder' 'settingsfile' = start main application and process
%                             files in 'folder' with 'settingsfile'
%
% For processing files without user interaction create a settings.mat-file
% (always created when processing single or multiple files), store that file
% in the 'Search Path' and start main application with parameter 'folder'
% and optionaly 'settings-file' (eg path to settings.mat)
%
% For automatic mode of the application (e.g. a folder is continously
% processed while users add new data) store a settings-file as
% 'settingsauto.mat' in same folder as application and restart

function TAPEEG(varargin)

if ismcc | isdeployed
    Path = ctfroot;
    tmp = strfind(Path,filesep);
    if ~isempty(tmp) & tmp(end) < length(Path)
        Path = Path(1:tmp(end));
    elseif ~isempty(tmp) & length(tmp) > 1
        Path = Path(1:tmp(end-1));
    else
        Path = pwd;
    end
else
    Path = mfilename('fullpath');
    tmp = strfind(Path,filesep);
    if ~isempty(tmp)
        Path = Path(1:tmp(end));
    else
        Path = pwd;
    end
    if exist(fullfile(Path,['java_libs' filesep 'poi-3.8-20120326.jar']),'file')
        javaaddpath(fullfile(Path,['java_libs' filesep 'poi-3.8-20120326.jar']));
        javaaddpath(fullfile(Path,['java_libs' filesep 'poi-ooxml-3.8-20120326.jar']));
        javaaddpath(fullfile(Path,['java_libs' filesep 'poi-ooxml-schemas-3.8-20120326.jar']));
        javaaddpath(fullfile(Path,['java_libs' filesep 'xmlbeans-2.3.0.jar']));
        javaaddpath(fullfile(Path,['java_libs' filesep 'dom4j-1.6.1.jar']));
    end
    if exist(fullfile(Path,['java_libs' filesep 'MFF-1.2.jar']),'file')
        javaaddpath(fullfile(Path,['java_libs' filesep 'MFF-1.2.jar']));
    end
    if exist('lab_read_data') ~= 2 %#ok<EXIST>
        addpath(Path);
        addpath(genpath(fullfile(Path,'functions')));
    end
end

global Main_Path VERSION
Main_Path = Path;
VERSION = 'v2.9';
clearvars Path

disp(['Start TAPEEG ' VERSION ' in ' Main_Path]);

if ~nargin
    % Set 'FixedWidthFontName' (only initial run)
    lab_selectfont_fixedwidth(false)
end

if ~nargin
    try %#ok<TRYNC>
        load(fullfile(Main_Path,'settingsauto.mat'));
    end
    if exist('cfg','var')
        cfg.SettingsFileAuto = true;
        cfg.settings_file = 'settingsauto.mat';
        cfg.settings_path = Main_Path;
        cfg.MAIN.log_path = fullfile(Main_Path,'TAPEEG.log');
        stayalive = 1; %#ok<NASGU>
    end
elseif nargin > 0
    searchfolder = varargin{1};
    if nargin > 1
        try
            load(varargin{2});
            if ~exist('cfg','var')
                disp('Abort, invalid settings-file')
                return
            end
        catch %#ok<CTCH>
            disp('Abort, settings-file not found')
            return
        end
    end
    if exist('cfg','var')
        if ~exist(searchfolder,'dir')
            disp('Abort, invalid search-path')
            return
        end
        if isfield(cfg,'SEARCH') & ~isempty(cfg.SEARCH) %#ok<NODEF>
            cfg.SEARCH.searchfolder = searchfolder;
            cfg.MAIN.log_path=fullfile(cfg.SEARCH.searchfolder,'TAPEEG.log');
            cfg.SettingsFileAuto = false;
            [~,~,~,cfg.settings_file] = lab_filename(varargin{2});
            cfg.settings_file = [cfg.settings_file '_run.mat'];
            cfg.settings_path = cfg.SEARCH.searchfolder;
        elseif isfield(cfg,'CollectConnect') & ~isempty(cfg.CollectConnect)
            cfg.CollectConnect.searchfolder = searchfolder;
        elseif isfield(cfg,'CollectGraph') & ~isempty(cfg.CollectGraph)
            cfg.CollectGraph.searchfolder = searchfolder;
        elseif isfield(cfg,'CollectFFT') & ~isempty(cfg.CollectFFT)
            cfg.CollectFFT.searchfolder = searchfolder;
        elseif isfield(cfg,'CollectMicro') & ~isempty(cfg.CollectMicro)
            cfg.CollectMicro.searchfolder = searchfolder;
        end
    else
        cfg.SEARCH.searchfolder = searchfolder;
        if exist(searchfolder,'dir')
            cfg.MAIN.log_path = fullfile(cfg.SEARCH.searchfolder,'TAPEEG.log');
            cfg.SettingsFileAuto = false;
            [~,~,~,cfg.settings_file] = lab_filename(varargin{2});
            cfg.settings_file = [cfg.settings_file '_run.mat'];
            cfg.settings_path = cfg.SEARCH.searchfolder;
        end
    end
end

if exist('stayalive','var') & exist('cfg','var') & isfield(cfg,'SEARCH') & isfield(cfg.SEARCH,'searchfolder') & exist(cfg.SEARCH.searchfolder,'dir')
    while 1 < 2
        try
            disp('')
            cfg = lab_tapeeg(cfg);
            disp('Wait 10 minutes and restart')
            pause(600);
        catch %#ok<CTCH>
            disp('Error processing files, try in 10 min again')
            pause(600);
        end
    end
elseif exist('stayalive','var') & exist('cfg','var') & isfield(cfg,'SEARCH')
    lab_tapeeg(cfg);
elseif exist('cfg','var') & isfield(cfg,'SEARCH') & isfield(cfg.SEARCH,'searchfolder') & exist(cfg.SEARCH.searchfolder,'dir')
    lab_tapeeg(cfg);
elseif exist('cfg','var') & isfield(cfg,'SEARCH') & isfield(cfg.SEARCH,'searchfolder') & strcmp(cfg.SEARCH.searchfolder,'tapeeg')
    lab_tapeeg;
elseif exist('cfg','var') & isfield(cfg,'SEARCH') & isfield(cfg.SEARCH,'searchfolder') & strcmp(cfg.SEARCH.searchfolder,'statistics')
    lab_main_statistics;
elseif exist('cfg','var') & isfield(cfg,'SEARCH') & isfield(cfg.SEARCH,'searchfolder') & strcmp(cfg.SEARCH.searchfolder,'plot')
    lab_main_plot;
elseif exist('cfg','var') & isfield(cfg,'CollectConnect') & isfield(cfg.CollectConnect,'searchfolder') & exist(cfg.CollectConnect.searchfolder,'dir')
    lab_collect_connectivity(cfg);
elseif exist('cfg','var') & isfield(cfg,'CollectGraph') & isfield(cfg.CollectGraph,'searchfolder') & exist(cfg.CollectGraph.searchfolder,'dir')
    lab_collect_graphanalysis(cfg);
elseif exist('cfg','var') & isfield(cfg,'CollectFFT') & isfield(cfg.CollectFFT,'searchfolder') & exist(cfg.CollectFFT.searchfolder,'dir')
    lab_collect_spectraldata(cfg);
elseif exist('cfg','var') & isfield(cfg,'CollectMicro') & isfield(cfg.CollectMicro,'searchfolder') & exist(cfg.CollectMicro.searchfolder,'dir')
    lab_collect_microstates(cfg);
else
    lab_show_start;
end