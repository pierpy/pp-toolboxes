% Search results of graphanalysis
%
% FilesAll = lab_collect_graphanalysis_search(cfg,skipselection)
%
% written by F. Hatz 2013

function FilesAll = lab_collect_graphanalysis_search(cfg,skipselection)

FilesAll = [];
if ~exist('skipselection','var')
    skipselection = false;
end
if ~exist('cfg','var')
    cfg = [];
end
if ~isfield(cfg,'CollectGraph') | ~isfield(cfg.CollectGraph,'searchfolder') | ~exist(cfg.CollectGraph.searchfolder,'dir')
    return
end

% Correct searchfolder if needed
searchfolder = cfg.CollectGraph.searchfolder;
if strcmp(searchfolder(end),filesep)
    searchfolder = searchfolder(1:end-1);
end
cd(searchfolder);

disp ('Search Files with results of graphanalysis')
searchstring = {'_GRAPH.mat'};
for i = 1:length(cfg.CollectGraph.includestring)
    searchstring{end+1,1} = ['+' cfg.CollectGraph.includestring{i}]; %#ok<AGROW>
end
for i = 1:length(cfg.CollectGraph.excludestring)
    searchstring{end+1,1} = ['|' cfg.CollectGraph.excludestring{i}]; %#ok<AGROW>
end
Files = lab_search(searchfolder,searchstring,skipselection);
if isempty(Files)
    return
end
Files = Files(:)';

Strings = {'PLI','dPLI','wPLI','wPLV','PLV','PLT', ...
    'AEC','EAC','SLc','SL',''};
for nstring = 1:size(Strings,2)
    if ~isempty(Strings{nstring})
        listtmp = Files(~cellfun(@isempty,strfind(upper(Files),upper(Strings{nstring}))));
    else
        listtmp = Files;
    end
    Files = setdiff(Files,listtmp);
    FilesAll(1,end+1).list = listtmp; %#ok<AGROW>
    if ~isempty(Strings{nstring})
        FilesAll(1,end).name = Strings{nstring};
    else
        FilesAll(1,end).name = 'GraphMeasure';
    end
    clearvars listbands listtmp tmp bands i
end

% select subjects
if skipselection == false
    Files = {};
    for i = 1:length(FilesAll)
        Files = cat(2,Files,FilesAll(i).list);
    end
    if isempty(Files)
        FilesAll = [];
        return
    end
    disp ('Select Files')
    selection = listdlg('PromptString','Files:','SelectionMode','multiple', ...
        'ListString',Files,'InitialValue',1:length(Files),'CancelString','None','ListSize',[450 400]);
    pause(0.2);
    if ~isempty(selection)
        Files = Files(1,selection);
        flags = zeros(1,length(FilesAll));
        for i = 1:length(FilesAll)
            Ftmp = intersect(Files,FilesAll(i).list);
            if ~isempty(Ftmp)
                FilesAll(i).list = Ftmp(:)'; %#ok<AGROW>
                flags(i) = 1;
            end
            clearvars Ftmp
        end
        if max(flags) == 1
            FilesAll = FilesAll(1,flags==1);
        else
            FilesAll = [];
        end
        clearvars selection
    else
        FilesAll = [];
    end
end

end