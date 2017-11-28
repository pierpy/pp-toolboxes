% Search and collect distances (from IS results)
%
% Files,skipprocessing] = lab_collect_distances_search(cfg,skipselection)
%
% written by F. Hatz 2014

function Files = lab_collect_distances_search(cfg,skipselection)

if ~exist('skipselection','var')
    skipselection = false;
end
if ~exist('cfg','var') | ~isfield(cfg,'CollectDistance') | ...
        ~isfield(cfg.CollectDistance,'searchfolder') | ...
        ~exist(cfg.CollectDistance.searchfolder,'dir')
    Files = {};
    return
end

% Correct searchfolder if needed
searchfolder = cfg.CollectDistance.searchfolder;
if strcmp(searchfolder(end),filesep)
    searchfolder = searchfolder(1:end-1);
end
cd(searchfolder);

disp('Search for result files (Distance.mat)')
Files = lab_search(searchfolder,'Distance.mat',skipselection);
Files = Files(:)';

% search for files to exclude
if ~isempty(Files)
    excludestring = cfg.CollectDistance.excludestring;
    for i = 1:size(excludestring,1)
        tmp = cellfun('isempty',strfind(Files,excludestring{i,1}));
        if ~isempty(tmp)
            Files = Files(1,tmp);
        else
            Files = [];
        end
    end
    clearvars excludestring i
end

% search for files to include
if ~isempty(Files)
    includestring = cfg.CollectDistance.includestring;
    for i = 1:size(includestring,1)
        tmp = ~cellfun('isempty',strfind(Files,includestring{i,1}));
        if ~isempty(tmp)
            Files = Files(1,tmp);
        else
            Files = [];
        end
    end
    clearvars includestring i
end

if ~isempty(Files) & skipselection == false
    selection = listdlg('PromptString','Select Files','SelectionMode','multiple', ...
        'ListString',Files,'InitialValue',1:length(Files),'CancelString','None','ListSize',[300 400]);
    if ~isempty(selection)
        Files = Files(selection);
        clearvars selection
    end
end


end