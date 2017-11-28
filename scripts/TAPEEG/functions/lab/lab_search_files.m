% Function to search files by defining searchfolder, search-strings and
% exclude-strings
%
% [calc,cfg] = lab_search_files(cfg,doextra)
%
% written by F. Hatz 2013

function [calc,cfg] = lab_search_files(cfg,doextra,dosort)

if ~exist('dosort','var')
    dosort = 1;
end
if ~exist('doextra','var')
    doextra = 0;
end
if ~exist('cfg','var')
    cfg = [];
end
if isfield(cfg,'MAIN') & isfield(cfg.MAIN,'auto') & cfg.MAIN.auto == 1
    skipselection = true;
else
    skipselection = false;
end

if ~isfield(cfg,'SEARCH') | ~isfield(cfg.SEARCH,'searchfolder')
    [cfg,skipprocessing] = lab_set_searchstrings(cfg,1,doextra);
    if skipprocessing == 1;
        calc.Filelist = [];
        return
    end
elseif ~isfield(cfg.SEARCH,'searchstring')
    [cfg,skipprocessing] = lab_set_searchstrings(cfg,0,doextra);
    if skipprocessing == 1;
        calc.Filelist = [];
        return
    end
end
if ~isfield(cfg.SEARCH,'excludestring')
    cfg.SEARCH.excludestring = {};
end
if ~isfield(cfg.SEARCH,'includestring')
    cfg.SEARCH.includestring = {};
end

% Correct searchfolder if needed
searchfolder = cfg.SEARCH.searchfolder;
if isempty(searchfolder) | ~ischar(searchfolder) | ~exist(searchfolder,'dir')
    calc.Filelist = [];
    return
end
if strcmp(searchfolder(end),filesep)
    searchfolder = searchfolder(1:end-1);
end

% Search for lists with processed files
if isfield(cfg,'settings_path') & ~isempty(cfg.settings_path)
    settings_path = cfg.settings_path;
else
    settings_path = searchfolder;
end
if exist(fullfile(settings_path,'ProcessedFiles.xls'),'file')
    if ispc
        [~,calc.Filelist_doneall] = xlsread(fullfile(settings_path,'ProcessedFiles.xls'));
    else
        [~,calc.Filelist_doneall] = xlsread(fullfile(settings_path,'ProcessedFiles.xls'),1,'','basic');
    end
    calc.Filelist_doneall= calc.Filelist_doneall(:,1)';
elseif exist(fullfile(settings_path,'ProcessedFiles.txt'),'file')
    calc.Filelist_doneall = importdata(fullfile(settings_path,'ProcessedFiles.txt'));
    calc.Filelist_doneall= calc.Filelist_doneall(:)';
else
    calc.Filelist_doneall = [];
end
if exist(fullfile(settings_path,'ProcessedFilesAll.xls'),'file')
    if ispc
        [~,calc.Filelist_done] = xlsread(fullfile(settings_path,'ProcessedFilesAll.xls'));
    else
        [~,calc.Filelist_done] = xlsread(fullfile(settings_path,'ProcessedFilesAll.xls'),1,'','basic');
    end
    calc.Filelist_done = calc.Filelist_done(:,1)';
elseif exist(fullfile(settings_path,'ProcessedFilesAll.txt'),'file')
    calc.Filelist_done = importdata(fullfile(settings_path,'ProcessedFilesAll.txt'));
    calc.Filelist_done= calc.Filelist_done(:)';
else
    calc.Filelist_done = [];
end

% correct for bad characters in folder names
List = dir(searchfolder);
if size(List,1) > 2
    for i = 3:size(List,1)
        if strcmp(List(i,1).name(end),'?')
            tmp = List(i,1).name(1:end-1);
            if ispc
                system(['move "' fullfile(searchfolder,[tmp '*']) '" "' fullfile(searchfolder,tmp) '"']);
            else
                system(['mv ''' fullfile(searchfolder,[tmp '*']) ''' ''' fullfile(searchfolder,tmp) '''']);
            end
        end
    end
end

% Search files
calc.Filelist = [];
skipICA = false;
doMFF = false;
searchstring = cat(1,cfg.SEARCH.searchstring(:),{'|.mrk';'|.vrb';'|part2of2_VE.fif';'|Spectra.mat';'|.info'});
if length(cfg.SEARCH.searchstring{1}) > 3 & strcmp(cfg.SEARCH.searchstring{1}(end-3:end),'.mrk')
    searchstring = setdiff(searchstring,{'|.mrk'},'stable');
end
for i = 1:length(cfg.SEARCH.searchstring)
    if strfind(searchstring{i},'.info')
        searchstring = setdiff(searchstring,{'|.info'},'stable');
    elseif strfind(searchstring{i},'.vrb')
        searchstring = setdiff(searchstring,{'|.vrb'},'stable');
    elseif strfind(searchstring{i},'.mff')
        searchstring{i} = 'signal1.bin';
        doMFF = true;
    elseif strfind(searchstring{i},'ICA.mat')
        skipICA = true;
    elseif strfind(searchstring{i},'.mat')
        searchstring = cat(1,searchstring,{'|ICA.mat'});
    elseif strcmp(searchstring{i},'.sef')
        searchstring = cat(1,searchstring,{'|ICA.sef';'|ICAtopo.sef';'|filt.sef';'|sef.vrb'});
    elseif strcmp(searchstring{i},'.edf')
        searchstring = cat(1,searchstring,{'|filt.edf';'|edf.vrb'});
    end
end
for i = 1:length(cfg.SEARCH.excludestring)
    searchstring{end+1,1} = ['|' cfg.SEARCH.excludestring{i}]; %#ok<AGROW>
end
for i = 1:length(cfg.SEARCH.includestring)
    searchstring{end+1,1} = ['+' cfg.SEARCH.includestring{i}]; %#ok<AGROW>
end
calc.Filelist = lab_search(searchfolder,searchstring,skipselection);
if doMFF == true
    Idx = find(~cellfun(@isempty,strfind(calc.Filelist,'signal1.bin')));
    for j = Idx
        calc.Filelist{j} = [calc.Filelist{j}(1:end-12) '.mff'];
    end
    clearvars j Idx
end

% Exclude processed files
if ~isempty(calc.Filelist_doneall)
    if isfield(cfg.SEARCH,'excludeprocessed') & cfg.SEARCH.excludeprocessed == true
        disp ('Exclude already processed files')
        calc.Filelist = setdiff(calc.Filelist,calc.Filelist_doneall);
        calc.Filelist = calc.Filelist(:)';
    end
end

% Include ICA results
if skipICA == false & isfield(cfg.SEARCH,'doICAresult') & cfg.SEARCH.doICAresult == true
    Filelist_ICA = lab_search(searchfolder,'exclude.txt',skipselection);
    for i = 1:length(Filelist_ICA)
        Filelist_ICA{i} = [Filelist_ICA{i}(1:end-12) '.mat'];
    end
    if isfield(cfg.SEARCH,'excludeprocessed') & cfg.SEARCH.excludeprocessed == true
        Filelist_ICA = setdiff(Filelist_ICA,calc.Filelist_doneall);
    end
    if ~isempty(Filelist_ICA)
        calc.Filelist = union(calc.Filelist,Filelist_ICA);
        calc.Filelist = calc.Filelist(:)';
    end
end

% Sort files by frequency bands
if dosort == 1
    disp ('Sort frequency bands')
    [calc.Filelist,calc.Freqbands] = lab_sort_frequencybands(calc.Filelist);
end

% Select files
if size(calc.Filelist,2) > 0 & skipselection == false
    disp ('Select Files')
    strlist = calc.Filelist;
    strdefault = 1:size(calc.Filelist,2);
    selection = listdlg('PromptString','Files:','SelectionMode','multiple', ...
        'ListString',strlist,'InitialValue',strdefault,'CancelString','None','ListSize',[450 400]);
    calc.Filelist = calc.Filelist(1,selection);
    if isfield(calc,'Freqbands')
        calc.Freqbands = calc.Freqbands(1,selection);
    end
    clearvars selection strlist strdefault
    pause(0.2);
end
