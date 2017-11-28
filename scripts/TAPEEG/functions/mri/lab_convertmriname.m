function MRI_file_out = lab_convertmriname(MRI_file)

[~,filepath,fileformat,filenameS] = lab_filename(MRI_file);
if strcmp(fileformat,'dcm') | strcmp(fileformat,'ima') | strcmp(fileformat,'img') | (isempty(fileformat) & ~isempty(filepath))
    tmp = strfind(filepath,filesep);
    if tmp(end) == length(filepath)
        filenameST = filepath(tmp(end-1)+1:tmp(end)-1);
        filepathT = filepath(1:tmp(end-1));
    else
        filenameST = filepath(tmp(end)+1:end);
        filepathT = filepath(1:tmp(end));
    end
    MRI_file_out = fullfile(filepathT,[filenameST '.hdr']);
    clearvars filenameST filepathT
else
    MRI_file_out = fullfile(filepath,[filenameS '.hdr']);
end
clearvars filepath filenameS fileformat

end