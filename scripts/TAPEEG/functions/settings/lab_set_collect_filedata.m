function [cfg,skipprocessing] = lab_set_collect_filedata(cfg)

disp ('Ask for collect file data settings')

skipprocessing = 0;

if ~exist('cfg','var') | ~isfield(cfg,'CollectFiles') | ~isfield(cfg.CollectFiles,'searchfolder')
    cfg.CollectFiles.searchfolder = pwd;
    cfg.CollectFiles.stringfilename = '';
    cfg.CollectFiles.stringpath = '';
    cfg.CollectFiles.estringfilename = '';
    cfg.CollectFiles.estringpath = '';
    cfg.CollectFiles.outputfolder = fullfile(pwd,'CollectFiles');
    cfg.CollectFiles.INTERPOLATE = [];
    cfg.CollectFiles.docopy = true;
    cfg.CollectFiles.doaverage = false;
    cfg.CollectFiles.grandaverage = false;
    cfg.CollectFiles.averagesearchstrings = '';
    cfg.CollectFiles.averagestring = '';
    cfg.CollectFiles.dogfp = false;
    cfg.CollectFiles.subjectname = 0;
end

if isfield(cfg,'EEG_file')
    Output = lab_prepare_subjectname(fullfile(cfg.EEG_filepath,cfg.EEG_file));
else
    Output = {[],[]};
end

Prompt = cell(0,2);
Formats = [];

Prompt(end+1,:) = {'Search Folder','searchfolder'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'dir';
Formats(end,1).size = 350;
Formats(end,1).callback = {@set_searchfolder,'@ALL','@ALL'};
Formats(end,1).span = [1 2];

Prompt{end+1,1} = '';
Formats(end+1,1).type = 'text';
Formats(end,1).span = [1 2];

Prompt{end+1,1} = 'Searchstrings Filename';
Formats(end+1,1).type = 'text';

Prompt(end+1,:) = {'','stringfilename'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'vector';
Formats(end,1).size = 300;

Prompt{end+1,1} = 'Searchstrings Filepath';
Formats(end+1,1).type = 'text';

Prompt(end+1,:) = {'','stringpath'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'vector';
Formats(end,1).size = 300;

Prompt{end+1,1} = 'Excludestrings Filename';
Formats(end+1,1).type = 'text';

Prompt(end+1,:) = {'','estringfilename'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'vector';
Formats(end,1).size = 300;

Prompt{end+1,1} = 'Excludestrings Filepath';
Formats(end+1,1).type = 'text';

Prompt(end+1,:) = {'','estringpath'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'vector';
Formats(end,1).size = 300;

Prompt{end+1,1} = '';
Formats(end+1,1).type = 'text';
Formats(end,1).span = [1 2];

Prompt(end+1,:) = {'Output folder','outputfolder'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'text';
Formats(end,1).size = 350;

Prompt{end+1,1} = '';
Formats(end+1,1).type = 'text';
Formats(end,1).span = [1 2];

Prompt(end+1,:) = {'Copy eeg/meg files','docopy'};
Formats(end+1,1).type = 'check';
Formats(end,1).span = [1 2];

Prompt{end+1,1} = '';
Formats(end+1,1).type = 'text';
Formats(end,1).span = [1 2];

Prompt(end+1,:) = {'Interpolate bad channels','INTERPOLATE'};
Formats(end+1,1).type = 'check';
Formats(end,1).span = [1 2];
Formats(end,1).callback = {@set_interpolate,'INTERPOLATE','INTERPOLATE'};

Prompt(end+1,:) = {'Average eeg/meg files','doaverage'};
Formats(end+1,1).type = 'check';
Formats(end,1).span = [1 2];

Prompt(end+1,:) = {'Grand Average','grandaverage'};
Formats(end+1,1).type = 'check';
Formats(end,1).span = [1 2];

Prompt{end+1,1} = 'Searchstrings for average';
Formats(end+1,1).type = 'text';

Prompt(end+1,:) = {'','averagesearchstrings'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'vector';
Formats(end,1).size = 300;

Prompt{end+1,1} = 'String for average files';
Formats(end+1,1).type = 'text';

Prompt(end+1,:) = {'','averagestring'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'vector';
Formats(end,1).size = 300;

Prompt{end+1,1} = '';
Formats(end+1,1).type = 'text';
Formats(end,1).span = [1 2];

Prompt(end+1,:) = {'Collect GFP','dogfp'};
Formats(end+1,1).type = 'check';
Formats(end,1).span = [1 2];

Prompt(end+1,:) = {'Number of underscores in subject name',''};
Formats(end+1,1).type = 'text';
Formats(end,1).span = [1 2];

if ~isempty(Output{1})
    Prompt(end+1,:) = {[Output{1} ' ' Output{2}],'subjectname'};
else
    Prompt(end+1,:) = {Output{2},'subjectname'};
end
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'integer';
Formats(end,1).limits = [-1 99];
Formats(end,1).size = 40;

[cfg.CollectFiles,Cancelled] = inputsdlg(Prompt,'Collect file data',Formats,cfg.CollectFiles);
if isempty(cfg.CollectFiles) | Cancelled == 1
    cfg.CollectFiles = [];
    skipprocessing = 1;
    return
else
    pause(0.2);
    if isfield(cfg.CollectFiles,'searchfolder') & exist(cfg.CollectFiles.searchfolder,'dir')
        save(fullfile(cfg.CollectFiles.searchfolder,'settings_files.mat'),'cfg','-v7.3');
    end
end

end

function settings = set_searchfolder(settings)
    if exist(fullfile(settings.searchfolder,'settings_files.mat'),'file')
        searchfolder = settings.searchfolder;
        load(fullfile(settings.searchfolder,'settings_files.mat'))
        if exist('cfg','var') & isfield(cfg,'CollectFiles') & ~isempty(cfg.CollectFiles)
            settings = cfg.CollectFiles;
            settings.searchfolder = searchfolder;
        end
    end
    if ~isempty(settings.searchfolder) & exist(settings.searchfolder,'dir') & isempty(settings.outputfolder)
        settings.outputfolder = fullfile(settings.searchfolder,'CollectFiles');
    else
        settings.outputfolder = '';
    end
end

function settings = set_interpolate(settings)
  Prompt = cell(0,2);
  Formats = [];
  
  Prompt(end+1,:) = {'Define bad channels','definebad'};
  Formats(end+1,1).type = 'edit';
  Formats(end,1).format = 'vector';
  Formats(end,1).limits = [0 inf];
  
  Prompt(end+1,:) = {'Detect bad','BAD'};
  Formats(end+1,1).type = 'edit';
  Formats(end,1).format = 'result';
  Formats(end,1).callback = {@lab_set_detect_bad,'BAD','BAD',[],[],0,0,0,0,1,0,1,1,0,1};
  Formats(end,1).size = [150 170];
  
  Prompt(end+1,:) = {'Interpolate bad channels','dointerpolate'};
  Formats(end+1,1).type = 'check';
  
  [settings,Cancelled] = inputsdlg(Prompt,'Interpolate bad channels',Formats,settings);
  if isempty(settings) | Cancelled == 1
      settings = [];
      return
  else
      pause(0.2);
  end
end