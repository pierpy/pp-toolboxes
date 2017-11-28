% Helper function for lab_collect_spectras
%
% [cfg,Filepaths_FFT,Filepaths_ISFFT] = lab_collect_spectraldata_search(cfg,skipselection)
%
% written by F. Hatz 2012

function Files = lab_collect_spectraldata_search(cfg,skipselection)

if ~exist('cfg','var')
    cfg = [];
end
if ~exist('skipselection','var')
    skipselection = false;
end
if ~isfield(cfg,'CollectFFT') | ~isfield(cfg.CollectFFT,'searchfolder') | ~exist(cfg.CollectFFT.searchfolder,'dir')
    Files = [];
    disp('   No valid folder selected')
    return
end

% Correct searchfolder if needed
searchfolder = cfg.CollectFFT.searchfolder;
if strcmp(searchfolder(end),filesep)
    searchfolder = searchfolder(1:end-1);
end
cd(searchfolder);

% Search files
disp('   Search for FFT-Results')
searchstring = {'Spectra.mat'};
for i = 1:length(cfg.CollectFFT.includestring)
    searchstring{end+1,1} = ['+' cfg.CollectFFT.includestring{i}]; %#ok<AGROW>
end
for i = 1:length(cfg.CollectFFT.excludestring)
    searchstring{end+1,1} = ['|' cfg.CollectFFT.excludestring{i}]; %#ok<AGROW>
end
if isfield(cfg.CollectFFT,'searchIS') & cfg.CollectFFT.searchIS == true
    searchstring{end+1,1} = '+_IS_';
end
Files = lab_search(searchfolder,searchstring,skipselection);
if isempty(Files)
    return
end
Files = Files(:)';

% Select files
if ~isempty(Files) & skipselection == false
    % Search for saved files in cfg
    disp ('Select Files')
    selection = listdlg('PromptString','Files:','SelectionMode','multiple', ...
        'ListString',Files,'InitialValue',1:length(Files),'CancelString','None','ListSize',[450 400]);
    if ~isempty(selection)
        Files = Files(1,selection);
    else
        Files = [];
    end
    clearvars selection strlist strdefault
    pause(0.2);
end
