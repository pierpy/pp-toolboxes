% Search Files in directory and subdirectory
%
% Files = lab_search(Searchpath,Pattern,Skipbar,Mode)
%
% Searchpath       full path
% Pattern          String or cell-arry with strings (wildcard * is accepted)
%                      String with +... = restrict results to ...
%                      String with |... = exclude results with ...
% Skipbar          (optional) 1 = disable progress bar
% Mode             1 = Search-Pattern is string with wildcards
%                  2 = Search-Pattern is regular expression
%
% Written by F. Hatz 2014

function [Files,Param] = lab_search(Searchpath,Pattern,Skipbar,NoSubFolder,Mode)
    
if ~exist('Mode','var')
    Mode = 1;
end
if ~exist('NoSubFolder','var')
    NoSubFolder = false;
end
if ~exist('Skipbar','var')
    Skipbar = false;
end
if ~exist('Searchpath','var')
    Searchpath = uigetdir('Select directory');
end
if ~exist('Pattern','var')
    Pattern{1} = '*';
elseif ischar(Pattern)
    tmp = Pattern;
    clearvars Pattern;
    Pattern = cellstr(tmp);
    clearvars tmp
end

% Prepare search pattern and extract include-pattern (+...) / exclude-patterns (|...)
NumPattern = length(Pattern);
Pexclude = [];
Pinclude = [];
for i = 1:NumPattern
    if length(Pattern{i}) > 1
        if strcmp(Pattern{i}(1),'|')
            Pexclude = [Pexclude i]; %#ok<AGROW>
        elseif strcmp(Pattern{i}(1),'+')
            Pinclude = [Pinclude i]; %#ok<AGROW>
        end
    end
end
if ~isempty(Pexclude)
    Exclude = Pattern(Pexclude);
    NumExclude = length(Exclude);
    for i = 1:NumExclude
        Exclude{i} = Exclude{i}(2:end);
    end
else
    NumExclude = 0;
end
if ~isempty(Pinclude)
    Include = Pattern(Pinclude);
    NumInclude = length(Include);
    for i = 1:NumInclude
        Include{i} = Include{i}(2:end);
    end
else
    NumInclude = 0;
end
Idx = setdiff(1:length(Pattern),union(Pinclude,Pexclude));
if ~isempty(Idx)
    Pattern = Pattern(Idx);
else
    Pattern = {'*'};
end
NumPattern = length(Pattern);
for i = 1:NumPattern
    Pattern{i} = regexprep(Pattern{i},{'/','\'},filesep);
    if Mode == 1
        Pattern{i} =  regexptranslate('wildcard',Pattern{i});
    end
end
clearvars Pexclude Pinclude i

% Find sub-directories
if NoSubFolder == false
    if Skipbar == false
        progressbar;
        progressbar('Search Files - lookup directory');
    end
    TmpList = genpath(Searchpath);
    Idx = strfind(TmpList,pathsep);
    NumI = length(Idx);
    DirList = cell(1,NumI+1);
    if isempty(Idx)
        DirList{1} = TmpList;
    else
        DirList{1} = TmpList(1:Idx(1)-1);
        for i = 1:NumI
            if i < NumI
                DirList{i+1} = TmpList(Idx(i)+1:Idx(i+1)-1);
            else
                DirList{i+1} = TmpList(Idx(i)+1:end);
            end
        end
    end
else
    DirList = cellstr(Searchpath);
end
NumDirs = length(DirList);

% Find Files in all directories
Files = {};
if Skipbar == false
    progressbar('Search Files');
end
Numdots = ceil(NumDirs/100);
for i = 1:NumDirs
    if mod(i,Numdots) == 0 & Skipbar == false
        progressbar(i/NumDirs);
    end
    Ftmp = dir(DirList{i});
    NumFiles = length(Ftmp);
    for j = 1:NumFiles
        if Ftmp(j).isdir ~= true
            if length(Ftmp(j).name) > 1 & strcmp(Ftmp(j).name(1:2),'._')
                delete(fullfile(DirList{i},Ftmp(j).name))
            elseif max(~cellfun(@isempty,regexpi(Ftmp(j).name,Pattern)))
                Files{1,end+1} = fullfile(DirList{i},Ftmp(j).name); %#ok<AGROW>
                if ~exist('Param','var')
                    Param = Ftmp(j);
                else
                    Param(1,end+1) = Ftmp(j); %#ok<AGROW>
                end
            end
        end
    end
end
if isempty(Files)
    Param = [];
    return
end

% Include Files
if NumInclude > 0
    if Skipbar == false
        progressbar('Include Files');
    end
    Idx = [];
    for i = 1:NumInclude
        Idx = union(Idx,find(~cellfun(@isempty,strfind(Files,Include{i}))));
    end
else
    Idx = 1:length(Files);
end

% Exclude Files
if NumExclude > 0
    if Skipbar == false
        progressbar('Exclude Files');
    end
    Idx2 = [];
    for i = 1:NumExclude
        Idx2 = union(Idx2,find(~cellfun(@isempty,strfind(Files,Exclude{i}))));
    end
    Idx = setdiff(Idx,Idx2);
end

if ~isempty(Idx)
    Files = Files(1,Idx);
    Param = Param(1,Idx);
else
    Files = {};
    Param = [];
end

if Skipbar == false
    progressbar(1);
end