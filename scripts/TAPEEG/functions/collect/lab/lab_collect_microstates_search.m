% Search files for lab_collect_microstates
%
% written by F. Hatz 2014

function Files = lab_collect_microstates_search(cfg,skipselection)

if ~exist('skipselection','var')
    skipselection = false;
end
if ~exist('cfg','var') | ~isfield(cfg,'CollectMicro') | ...
        ~isfield(cfg.CollectMicro,'searchfolder') | ...
        ~exist(cfg.CollectMicro.searchfolder,'dir')
    disp('   No valid folder selected')
    Files = [];
    return
end

% Correct searchfolder if needed
searchfolder = cfg.CollectMicro.searchfolder;
if strcmp(searchfolder(end),filesep)
    searchfolder = searchfolder(1:end-1);
end
cd(searchfolder);

% Search files
disp('Search for result files')
searchstring = {'Microstates.mat'};
for i = 1:length(cfg.CollectMicro.includestring)
    searchstring{end+1,1} = ['+' cfg.CollectMicro.includestring{i}]; %#ok<AGROW>
end
for i = 1:length(cfg.CollectMicro.excludestring)
    searchstring{end+1,1} = ['|' cfg.CollectMicro.excludestring{i}]; %#ok<AGROW>
end
Files = lab_search(searchfolder,searchstring,skipselection);
if isempty(Files)
    return
end
Files = Files(:)';

% Select files
if skipselection == false
    if isempty(Files)
        Files = [];
        return
    end
    disp ('Select Files')
    strlist = Files;
    strdefault = 1:size(Files,2);
    selection = listdlg('PromptString','Files:','SelectionMode','multiple', ...
        'ListString',strlist,'InitialValue',strdefault,'CancelString','None','ListSize',[450 400]);
    pause(0.2);
    if ~isempty(selection)
        Files = Files(1,selection);
    else
        Files = [];
        return
    end
    clearvars selection strlist strdefault
end
