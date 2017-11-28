% Search files for lab_collect_connectivity
%
% written by F. Hatz 2012

function [cfg,Files,skipprocessing] = lab_collect_connectivity_search_random(cfg)

skipprocessing = 0;

if ~exist('cfg','var')
    cfg = [];
end

if ~isfield(cfg,'CollectConnect') | ~isfield(cfg.CollectConnect,'searchfolder') | ~exist(cfg.CollectConnect.searchfolder,'dir')
    [cfg,skipprocessing] = lab_set_collect_connectivity(cfg);
    if skipprocessing == 1
        Files = [];
        return
    end
end

% Correct searchfolder if needed
searchfolder = cfg.CollectConnect.searchfolder;
if strcmp(searchfolder(end),filesep)
    searchfolder = searchfolder(1:end-1);
end
cd(searchfolder);

% Search files
searchstring = {'RandPhase_Conn_F*.mat'};
for i = 1:length(cfg.CollectConnect.includestring)
    searchstring{end+1,1} = ['+' cfg.CollectConnect.includestring{i}]; %#ok<AGROW>
end
for i = 1:length(cfg.CollectConnect.excludestring)
    searchstring{end+1,1} = ['|' cfg.CollectConnect.excludestring{i}]; %#ok<AGROW>
end
Filelist = lab_search(searchfolder,searchstring,true);
if isempty(Filelist)
    return
end
Filelist = Filelist(:)';

% Find files calculated for random phases
FilelistRand = [];

% Select files
disp ('Select Files')
strlist = Filelist;
strdefault = 1:size(Filelist,2);
selection = listdlg('PromptString','Files:','SelectionMode','multiple', ...
    'ListString',strlist,'InitialValue',strdefault,'CancelString','None','ListSize',[450 400]);
pause(0.2);
if ~isempty(selection)
    Filelist = Filelist(1,selection);
    if ~isempty(FilelistRand)
        FilelistRand = FilelistRand(1,selection);
    end
else
    Files = [];
    return
end
clearvars selection strlist strdefault

bands = [];
for i = 1:size(Filelist,2)
    tmp = strfind(Filelist{1,i},'_');
    tmp = Filelist{1,i}(tmp(end-1)+1:end-4);
    if ~strcmp(bands,tmp)
        bands = [bands cellstr(tmp)];
    end
    listbands(i) = find(strcmp(bands,tmp)); 
end

j = 1;
for i = 1:size(bands,2)
    Files(1,j).list = Filelist(1,listbands==i);
    if ~isempty(FilelistRand)
        Files(1,j).listrand = FilelistRand(1,listbands==i);
    end
    Files(1,j).name = ['Conn_' bands{1,i}];
    j = j+1;
end
clearvars listbands Filelist tmp bands i
