% Script to create for every eeg/meg-file a folder using patient-name and
% move the file inside
%
% written by F. Hatz 2012

function lab_change_filenames

[calc,cfg] = lab_search_files;
Filelist = calc.Filelist;
if isempty(Filelist)
    return
end

[Output] = lab_prepare_subjectname(Filelist{1});
tmp = strfind(Output{1},'^-');
tmp = setdiff(1:length(Output{1}),union(tmp,tmp+1));
Output{1} = Output{1}(tmp);
cfg.subjectname1 = 1;
cfg.subjectname2 = 0;

Prompt = cell(0,2);
Formats = [];

if ~isempty(Output{1})
    Prompt(end+1,:) = {'Number of folders to include',''};
    Formats(end+1,1).type = 'text';
    Prompt(end+1,:) = {Output{1},'subjectname1'};
    Formats(end+1,1).type = 'edit';
    Formats(end,1).format = 'integer';
    Formats(end,1).limits = [1 99];
    Formats(end,1).size = 40;
end

Prompt(end+1,:) = {'Number of underscores in subject name',''};
Formats(end+1,1).type = 'text';
Prompt(end+1,:) = {Output{2},'subjectname2'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'integer';
Formats(end,1).limits = [0 99];
Formats(end,1).size = 40;

[cfg,Cancelled] = inputsdlg(Prompt,'Subject name',Formats,cfg);
if Cancelled == 1
    return
else
    pause(0.2);
end


for i=1:length(Filelist)
    [Filename,Filepath,Format] = lab_filename(Filelist{i});
    if ~isempty(strfind(Filename,'_matrix.txt'))
        Format = '_matrix.txt';
    else
        Format = ['.' Format]; %#ok<AGROW>
    end
    cfg.subjectname = cfg.subjectname2;
    Filename = lab_subjectname(Filelist{i},cfg);
    for j = 1:cfg.subjectname1
        cfg.subjectname = -j;
        Filename = [lab_subjectname(Filelist{i},cfg) '_' Filename]; %#ok<AGROW>
    end
    movefile(Filelist{i},fullfile(Filepath,[Filename Format]));
end

end