% Write Covariance file .cov
%
% lab_write_cov(filename,header)
%
% written by F. Hatz 2014

function lab_write_cov(Filename,header)

if isnumeric(header) & size(header,1) == size(header,2)
    COV = header;
elseif isstruct(header) & isfield(header,'cov') & ~isempty(header.cov)
    COV = header.cov;
else
    return
end

[~,Filepath,~,FilenameS] = lab_filename(Filename);
dlmwrite(fullfile(Filepath,[FilenameS '.cov']),COV,'delimiter','\t','precision', 6);

end
