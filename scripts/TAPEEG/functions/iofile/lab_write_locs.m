function filename = lab_write_locs(filename,locs,badchans)

if isfield(locs,'locs')
    header = locs;
    locs = locs.locs;
end
if ~isnumeric(locs) & ~isfield(locs,'x') & ~isfield(locs,'chanpos') & ~isfield(locs,'coilpos')
    filename = [];
    return
end


[filename,filepath,format,filenameS] = lab_filename(filename);

if strcmp(format,'els') & exist('badchans','var')
    if exist('header','var')
        lab_write_els(fullfile(filepath,filename),header,badchans);
    else
        lab_write_els(fullfile(filepath,filename),locs,badchans);
    end
elseif exist(['lab_write_' format]) == 2 %#ok<EXIST>
    eval(['lab_write_' format '(fullfile(filepath,filename),locs);']);
else
    disp(['    Error writing loc-file, ' format ' not supported'])
end

if isfield(locs,'grad')
    grad = locs.grad; %#ok<NASGU>
    if isfield(locs,'digits')
        digits = locs.digits; %#ok<NASGU>
        save(fullfile(filepath,[filenameS '.grad']),'grad','digits');
    else
        save(fullfile(filepath,[filenameS '.grad']),'grad');
    end
end

