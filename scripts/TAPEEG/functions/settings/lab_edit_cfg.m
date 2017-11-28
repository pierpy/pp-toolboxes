% Show dialg for settings of TAPEEG
%
% written by F. Hatz 2012

function [cfg,header,skipprocessing] = lab_edit_cfg(cfg,header,nosearch,icaback,skipfilt,Title)

global THeader TData
skipprocessing = 0;

Filename = '';
if ~exist('Title','var') | isempty(Title)
    Title = 'Edit settings';
end
if ~exist('skipfilt','var')
    skipfilt = 0;
end
if ~exist('icaback','var')
    icaback = 0;
end
if ~exist('nosearch','var')
    nosearch = 0;
end
if ~exist('header','var')
    header = [];
elseif ischar(header)
    if exist(header,'file')
        Filename = header;
    end
    header = [];
end
if ~exist('cfg','var')
    cfg = [];
end

TData = [];
THeader = header;
if isempty(THeader)
    if isempty(Filename) & isfield(cfg,'SEARCH') & isfield(cfg.SEARCH,'searchfolder') & ...
            isfield(cfg.SEARCH,'searchstring')
        cfgtmp = cfg;
        cfgtmp.MAIN.auto = 1;
        calc = lab_search_files(cfgtmp);
        if ~isempty(calc.Filelist)
            Filename = calc.Filelist{1,1};
        end
        clearvars calc
    end
    if isempty(Filename)
        [Filename,Filepath] = uigetfile('*.*','Select file to preload header');
        if Filename ~= 0
            Filename = fullfile(Filepath,Filename);
        end
        clearvars Filepath
    end
    if ~isempty(Filename) & ischar(Filename) & exist(Filename,'file')
        if strcmp(Filename(end-3:end),'.mat')
            THeader = load_ICAmat(Filename);
            icaback = 1;
        else
            disp(['    preload header (' Filename ')'])
            [TData,THeader,cfg] = lab_read_data(Filename,cfg,true,[],true);
            if ~isfield(cfg,'SEARCH') | ~isfield(cfg.SEARCH,'searchfolder') | isempty(cfg.SEARCH.searchfolder)
                [~,cfg.SEARCH.searchfolder] = lab_filename(Filename);
            end
            if ~isfield(cfg,'SEARCH') | ~isfield(cfg.SEARCH,'searchstring') | isempty(cfg.SEARCH.searchstring)
                [~,~,cfg.SEARCH.searchstring{1}] = lab_filename(cfg.EEG_file);
                cfg.SEARCH.searchstring{1} = ['.' cfg.SEARCH.searchstring{1}];
            end
        end
    end
end

if isfield(cfg,'EEG_file')
    Output = lab_prepare_subjectname(fullfile(cfg.EEG_filepath,cfg.EEG_file));
else
    Output = {[],[]};
end
if ~isfield(cfg,'subjectname') | isempty(cfg.subjectname)
    if ~isempty(Output{1})
        cfg.subjectname = -1;
    else
        cfg.subjectname = 0;
    end
end
if ~isfield(cfg,'SEG') | ~isfield(cfg.SEG,'select')
    cfg.SEG.select = 'Select complete file';
end
if isfield(cfg,'SEG') & isfield(cfg.SEG,'MARK') & ~isempty(cfg.SEG.MARK) & ...
        (~isfield(cfg,'MARK') | isempty(cfg.MARK))
    cfg.MARK = cfg.SEG.MARK;
end
if isfield(THeader,'locs') & ~isempty(THeader.locs)
    if ~isfield(cfg,'LOCS') | ~isfield(cfg.LOCS,'locs') | isempty(cfg.LOCS.locs)
        cfg.LOCS.filelocs = true;
    end
else
    cfg.LOCS.filelocs = false;
    if ~isfield(cfg.LOCS,'locs') | isempty(cfg.LOCS.locs)
        cfg.LOCS = [];
    end
end

if isfield(cfg,'EXTRA') & isfield(cfg.EXTRA,'numdatachans') & ...
        isfield(THeader,'numdatachannels') & THeader.numdatachannels ~= cfg.EXTRA.numdatachans
    cfg.EXTRA = [];
    cfg.exclude = [];
end
if ~isfield(cfg,'EXTRA') | ~isfield(cfg.EXTRA,'numdatachans')
    cfg = lab_set_extrachannels(cfg,THeader,TData,1);
end

Prompt = cell(0,2);
Formats = [];

if nosearch == 0
    if ~isfield(cfg,'SEARCH') | ~isfield(cfg.SEARCH,'searchfolder')
        cfg.SEARCH.searchfolder = pwd;
    end
    if ~isfield(cfg.SEARCH,'searchstring')
        cfg.SEARCH.searchstring = {'.sef'};
    end
    if ~isfield(cfg.SEARCH,'searchmode')
        cfg.SEARCH.searchmode = 'Cartool (.sef)';
    end
    if ~isfield(cfg.SEARCH,'excludestring')
        cfg.SEARCH.excludestring = {''};
    end
    Prompt(end+1,:) = {'Search Files','SEARCH'};
    Formats(end+1,1).type = 'check';
    Formats(end,1).callback = {@lab_set_searchstrings,'@ALL','@ALL',1,1};
    
    Formats(end+1,1).type = 'none';
    Formats(end,1).span = [1 4];
end

Prompt(end+1,:) = {'Number of underscores in subject name',''};
Formats(end+1,1).type = 'text';
Formats(end,1).span = [1 4];

if ~isempty(Output{1})
    Prompt(end+1,:) = {[Output{1} ' ' Output{2}],'subjectname'};
else
    Prompt(end+1,:) = {Output{2},'subjectname'};
end
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'integer';
Formats(end,1).limits = [-99 99];
Formats(end,1).size = 40;

Formats(end+1,1).type = 'none';
Formats(end,1).span = [1 4];

if icaback ~= 1
    Prompt(end+1,:) = {'File reading',[]};
    Formats(end+1,1).type = 'text';
    Formats(end,1).size = [-1 0];
    Formats(end,1).span = [1 4]; % item is 1 field x 4 fields
    
    Prompt(end+1,:) = {'Segment reading','FILEREADING'};
    Formats(end+1,1).type = 'check';
    Formats(end,1).callback = {@lab_set_define_filereading,'@ALL','@ALL'};

    Prompt(end+1,:) = {'Locs','LOCS'};
    Formats(end+1,1).type = 'check';
    Formats(end,1).callback = {@set_locs,'@ALL','@ALL'};
    
    Prompt(end+1,:) = {'Extra channels','EXTRA'};
    Formats(end+1,1).type = 'check';
    Formats(end,1).callback = {@set_extrachannels,'@ALL','@ALL'};
    
    Prompt(end+1,:) = {'Add Ref-Channel','ADDREF'};
    Formats(end+1,1).type = 'check';
    Formats(end,1).callback = {@set_addref,'@ALL','@ALL'};
end

Formats(end+1,1).type = 'none';
Formats(end,1).span = [1 4];

Prompt(end+1,:) = {'Settings',[]};
Formats(end+1,1).type = 'text';
Formats(end,1).size = [-1 0];
Formats(end,1).span = [1 4]; % item is 1 field x 4 fields

Prompt(end+1,:) = {'Edit Markers','MARK'};
Formats(end+1,1).type = 'check';
Formats(end,1).callback = {@lab_set_edit_markers,'@ALL','@ALL'};

Prompt(end+1,:) = {'Create segments','SEG'};
Formats(end+1,1).type = 'check';
Formats(end,1).callback = {@lab_set_define_segments,'@ALL','@ALL'};

Prompt(end+1,:) = {'Exclude channels','exclude'};
Formats(end+1,1).type = 'check';
Formats(end,1).callback = {@lab_set_exclude,'@ALL','@ALL'};

Formats(end+1,1).type = 'none';
Formats(end,1).span = [1 4];

if icaback == 1
    [cfg,Prompt,Formats] = lab_set_preprocessingICA(cfg,THeader,Prompt,Formats);
else
    [cfg,Prompt,Formats] = lab_set_preprocessing(cfg,THeader,skipfilt,0,Prompt,Formats);
    if isfield(cfg,'ICA') & ~isempty(cfg.ICA) & isfield(cfg,'ICABACK') & ~isempty(cfg.ICABACK)
        cfg.ICA.ICABACK = cfg.ICABACK;
    end
end

Formats(end+1,1).type = 'none';
Formats(end,1).span = [1 4];

[cfg,Prompt,Formats] = lab_set_postprocessing(cfg,THeader,TData,Prompt,Formats);

[cfg,Cancelled] = inputsdlg(Prompt,Title,Formats,cfg);
if isempty(cfg) | Cancelled == 1
    cfg = [];
    skipprocessing = 1;
elseif isfield(cfg,'ICA') & ~isempty(cfg.ICA) & isfield(cfg.ICA,'ICABACK') & ~isempty(cfg.ICA.ICABACK)
    cfg.ICABACK = cfg.ICA.ICABACK;
end
pause(0.2);
header = THeader;
THeader = [];
TData = [];

    function settings = set_locs(settings)
        if isfield(settings.LOCS,'matchlocs') & settings.LOCS.matchlocs == true
            domatch = true;
        else
            domatch = false;
        end
        settings = lab_set_locs(settings);
        if isfield(settings.LOCS,'matchlocs') & settings.LOCS.matchlocs == true
            domatch2 = true;
        else
            domatch2 = false;
        end
        if ~isempty(Filename) & ischar(Filename) & exist(Filename,'file')
            if domatch ~= domatch2
                settings.EXTRA = [];
                [TData,THeader,settings] = lab_read_data(Filename,settings,true);
                settings = lab_set_extrachannels(settings,THeader,TData,1);
            else
                [TData,THeader,settings] = lab_read_data(Filename,settings,true);
            end
        elseif isfield(settings,'LOCS') & isfield(settings.LOCS,'locs')
            THeader.locs = settings.LOCS.locs;
        end
    end
    
    function settings = set_extrachannels(settings)
        if isfield(settings.EXTRA,'auxmethod') & isfield(settings.EXTRA,'numdatachans')
            Auxmethod = settings.EXTRA.auxmethod;
            Ref_chan = settings.EXTRA.ref_chan;
            Numdatachans = settings.EXTRA.numdatachans;
        else
            Auxmethod = '';
            Ref_chan = -1;
            Numdatachans =-1;
        end
        if isfield(settings.EXTRA,'reduceAAL')
            ReduceAAL = settings.EXTRA.reduceAAL;
        else
            ReduceAAL = -1;
        end
        settings = lab_set_extrachannels(settings,THeader,TData);
        if isfield(settings.EXTRA,'auxmethod') & isfield(settings.EXTRA,'numdatachans')
            Auxmethod2 = settings.EXTRA.auxmethod;
            Ref_chan2 = settings.EXTRA.ref_chan;
            Numdatachans2 = settings.EXTRA.numdatachans;
        else
            Auxmethod2 = '';
            Ref_chan2 = -1;
            Numdatachans2 = -1;
        end
        if isfield(settings.EXTRA,'reduceAAL')
            ReduceAAL2 = settings.EXTRA.reduceAAL;
        else
            ReduceAAL2 = -1;
        end
        if ~isempty(Filename) & ischar(Filename) & exist(Filename,'file') & ...
                (~strcmp(Auxmethod,Auxmethod2) | Ref_chan ~= Ref_chan2 | ...
                Numdatachans ~= Numdatachans2 | ReduceAAL ~= ReduceAAL2)
            [TData,THeader,settings] = lab_read_data(Filename,settings,true,[],true);
        end
    end
    
    function settings = set_addref(settings)
        if isfield(settings.ADDREF,'name') & ~isempty(settings.ADDREF.name)
            doref = settings.ADDREF.name;
        else
            doref = '';
        end
        settings = lab_set_addref(settings);
        if isfield(settings.ADDREF,'name') & ~isempty(settings.ADDREF.name)
            doref2 = settings.ADDREF.name;
        else
            doref2 = '';
        end
        if ~isempty(Filename) & ischar(Filename) & exist(Filename,'file') & ~strcmp(doref,doref2)
            settings.EXTRA = [];
            [TData,THeader,settings] = lab_read_data(Filename,settings,true,[],true);
            settings = lab_set_extrachannels(settings,THeader,TData,1);
        end
    end

end

function header = load_ICAmat(Filename)

try %#ok<TRYNC>
    load(Filename);
end
if ~exist('header','var')
    header = [];
end

end
