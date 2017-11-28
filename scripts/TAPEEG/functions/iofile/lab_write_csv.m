function lab_write_csv(Filename,Table)
    
[~,Filepath,~,FilenameS] = lab_filename(Filename);
header.main = {};
if min(cellfun(@isnumeric,Table(:,1))) == 0 & min(cellfun(@isnumeric,Table(1,:))) == 0
    header.subjects = Table(1,2:end);
    header.measures = Table(2:end,1)';
    header.main = Table(1,1);
    Table = Table(2:end,2:end);
elseif min(cellfun(@isnumeric,Table(1,:))) == 0
    header.subjects = Table(1,:);
    header.measures = {};
    Table = Table(2:end,:);
elseif min(cellfun(@isnumeric,Table(:,1))) == 0
    header.measures = Table(:,1)';
    header.subjects = {};
    Table = Table(:,2:end);
end
tmp = cellfun(@isnumeric,Table);
Nx = min(tmp,[],2) == 1;
Ny = min(tmp,[],1) == 1;
Table = cell2mat(Table(Nx,Ny));
if ~isempty(header.measures)
    header.measures = header.measures(Nx);
end
if ~isempty(header.subjects)
    header.subjects = header.subjects(Ny); %#ok<STRNU>
end
FileOut = fullfile(Filepath,[FilenameS '.csv']);
if exist(FileOut,'file')
    delete(FileOut);
end
csvwrite(FileOut,Table);

FileOut = fullfile(Filepath,[FilenameS '.head']);
if exist(FileOut,'file')
    delete(FileOut);
end
save(FileOut,'header');

end