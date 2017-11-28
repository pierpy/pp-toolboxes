% Script to create for every eeg/meg-file a folder using patient-name and
% move the file inside
%
% written by F. Hatz 2012

function lab_create_patientfolders

[calc,cfg] = lab_search_files;
Filelist = calc.Filelist;

cfg.newdir = false;
cfg.directory = '';

Prompt = cell(0,2);
Formats = [];

Prompt(end+1,:) = {'New folder for files' 'newdir'};
Formats(end+1,1).type = 'check';
Formats(end+1,1).enable = 'inactive';

Prompt(end+1,:) = {'Folder' 'directory'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'dir';
Formats(end,1).size = [-1 0];
Formats(end,1).callback = {@set_directory,'@ALL','@ALL'};

[cfg,Cancelled] = inputsdlg(Prompt,'Create patient folders',Formats,cfg);
if Cancelled == 1 | isempty(cfg)
    return
else
    pause(0.2);
end

for i=1:length(Filelist)
    [Filename,Filepath,~,FilenameS] = lab_filename(Filelist{i});
    [subjectname,cfg,skipprocessing] = lab_subjectname(Filename,cfg);
    if skipprocessing == 1
        return
    else
        pause(0.2);
    end
    
    if cfg.newdir == true & exist(cfg.directory,'dir')
        FilepathNew = fullfile(cfg.directory,subjectname);
    else
        FilepathNew = fullfile(Filepath,subjectname);
    end
    warning off %#ok<WNOFF>
    mkdir(FilepathNew);
    warning off %#ok<WNOFF>
    if strcmp(Filename(end-2:end),'mff') & exist(fullfile(Filepath,FilenameS),'dir')
        movefile(fullfile(Filepath,FilenameS),fullfile(FilepathNew,FilenameS));
    else
        movefile(fullfile(Filepath,Filename),fullfile(FilepathNew,Filename));
    end
    if exist(fullfile(Filepath,[Filename '.mrk']),'file')
        movefile(fullfile(Filepath,[Filename '.mrk']),fullfile(FilepathNew,[Filename '.mrk']));
    end
    if exist(fullfile(Filepath,[FilenameS '.mrk']),'file')
        movefile(fullfile(Filepath,[FilenameS '.mrk']),fullfile(FilepathNew,[FilenameS '.mrk']));
    end
    if exist(fullfile(Filepath,[FilenameS '.info']),'file')
        movefile(fullfile(Filepath,[FilenameS '.info']),fullfile(FilepathNew,[FilenameS '.info']));
    end
    if exist(fullfile(Filepath,[FilenameS '.els']),'file')
        movefile(fullfile(Filepath,[FilenameS '.els']),fullfile(FilepathNew,[FilenameS '.els']));
    end
    if exist(fullfile(Filepath,[FilenameS '.grad']),'file')
        movefile(fullfile(Filepath,[FilenameS '.grad']),fullfile(FilepathNew,[FilenameS '.grad']));
    end
end

end

function settings = set_directory(settings)

if ~isempty(settings.directory)
    settings.newdir = true;
else
    settings.newdir = false;
end

end