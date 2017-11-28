% Search files for lab_collect_connectivity
%
% written by F. Hatz 2012

function Files = lab_collect_connectivity_search(cfg,skipselection)

Files = [];
if ~exist('skipselection','var')
    skipselection = false;
end
if ~exist('cfg','var') | ~isfield(cfg,'CollectConnect') | ~isfield(cfg.CollectConnect,'searchfolder') | ~exist(cfg.CollectConnect.searchfolder,'dir')
    disp('No files with Connectivity-Results found - invalid search path')
    Files = [];
    return
end

% Correct searchfolder if needed
searchfolder = cfg.CollectConnect.searchfolder;
if strcmp(searchfolder(end),filesep)
    searchfolder = searchfolder(1:end-1);
end
cd(searchfolder);

% Search files
disp('Search for connectivity result files')
searchstring = {'Conn_F*.mat'};
for i = 1:length(cfg.CollectConnect.includestring)
    searchstring{end+1,1} = ['+' cfg.CollectConnect.includestring{i}]; %#ok<AGROW>
end
for i = 1:length(cfg.CollectConnect.excludestring)
    searchstring{end+1,1} = ['|' cfg.CollectConnect.excludestring{i}]; %#ok<AGROW>
end
Filelist = lab_search(searchfolder,searchstring,skipselection);
if isempty(Filelist)
    return
end
[Filelist,Freqbands] = lab_sort_frequencybands(Filelist);
Filelist = Filelist(:)';

% Find files calculated for random phases
Filenames = Filelist;
for i = 1:length(Filelist)
    Filenames{1,i} = lab_filename(Filelist{1,i});
end
T = strfind(Filenames,'RandPhase');
tmp = cellfun('isempty',T);
tmp2 = ~cellfun('isempty',T);
if ~isempty(tmp)
    Filelisttmp = Filelist;
    Filelist = Filelist(1,tmp);
    if ~isempty(tmp2)
        Filelisttmp = Filelisttmp(1,tmp2);
        FilelistRand = cell(1,length(Filelist));
        for i = 1:length(Filelisttmp)
            N = find(strcmp(Filelist(1,:),regexprep(Filelisttmp{1,i},'RandPhase_','')));
            if ~isempty(N)
                FilelistRand(1,N) = Filelisttmp(1,i);
            end
        end
        clearvars i N
    else
        FilelistRand = [];
    end
    clearvars Filelisttmp
else
    Files = [];
    return
end
clearvars T tmp tmp2

% Select files
if skipselection == false
    if isempty(Filelist)
        Files = [];
        return
    end
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
        Freqbands = Freqbands(1,selection);
    else
        Files = [];
        return
    end
    clearvars selection strlist strdefault
end

[Freqband,~,Fsort] = unique(Freqbands,'stable');
Files = [];
for i = 1:length(Freqband)
    if ~isempty(Fsort==i)
        Files(1,end+1).list = Filelist(1,Fsort==i); %#ok<AGROW>
        if ~isempty(FilelistRand)
            Files(1,end).listrand = FilelistRand(1,Fsort==i);
        end
        Files(1,end).name = ['Conn_' Freqband{1,i}];
        Files(1,end).freqband = Freqband{1,i};
    end
end
clearvars listbands Filelist tmp bands i
