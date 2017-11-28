function lab_write_matrix(Filename,matrix)

[~,Filepath,~,FilenameS] = lab_filename(Filename);

dlmwrite(fullfile(Filepath,[FilenameS '.txt']),matrix,'delimiter','\t','precision', 6);