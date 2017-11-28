function lab_write_xls(Filename,xlsout,Sheet)

global doxls

[~,Filepath,Format,FilenameS] = lab_filename(Filename);
if isempty(Format) | ~strcmp(Format,'xls') | ~strcmp(Format,'xlsx')
    Format = 'xlsx';
end
if length(FilenameS) > 80
    FilenameS = FilenameS(1:80);
end
if isempty(doxls)
    doxls = false;
    try %#ok<TRYNC>
        e = actxserver ('Excel.Application');
        if ~isempty(e) & ishandle(e)
            doxls = true;
        end
        e.Quit
        e.delete
        clear e
    end
end

if isnumeric(xlsout)
    xlsout = num2cell(xlsout);
end

if ~exist('Sheet','var') | ~ischar(Sheet)
    Sheet = '';
end

if doxls == true
    if length(xlsout(:)) > 2500000
        disp('Write csv-file, data array to large (>2''500''000 items)');
        FileOut = fullfile(Filepath,[FilenameS '.csv']);
        if exist(FileOut,'file')
            delete(FileOut);
        end
        lab_write_csv(FileOut,xlsout);
    elseif size(xlsout,2) > 5000
        disp('Write csv-file, to many columns (>5''000)');
        FileOut = fullfile(Filepath,[FilenameS '.csv']);
        if exist(FileOut,'file')
            delete(FileOut);
        end
        lab_write_csv(FileOut,xlsout);
    elseif length(xlsout(:)) > 100000 | ~exist('xlwrite') %#ok<EXIST>
        FileOut = fullfile(Filepath,[FilenameS '.' Format]);
        if exist(FileOut,'file') & isempty(Sheet)
            delete(FileOut);
        end
        if ~isempty(Sheet)
            xlswrite(FileOut,xlsout,Sheet);
        else
            xlswrite(FileOut,xlsout);
        end
    else
        FileOut = fullfile(Filepath,[FilenameS '.' Format]);
        if exist(FileOut,'file') & isempty(Sheet)
            delete(FileOut);
        end
        if ~isempty(Sheet)
            xlwrite(FileOut,xlsout,Sheet);
        else
            xlwrite(FileOut,xlsout);
        end
    end
else
    if length(xlsout(:)) > 100000 | ~exist('xlwrite') %#ok<EXIST>
        disp('Write csv-file, data array to large (>100''000 items)');
        FileOut = fullfile(Filepath,[FilenameS '.csv']);
        if exist(FileOut,'file')
            delete(FileOut);
        end
        lab_write_csv(FileOut,xlsout);
    else
        FileOut = fullfile(Filepath,[FilenameS '.' Format]);
        if exist(FileOut,'file') & isempty(Sheet)
            delete(FileOut);
        end
        if ~isempty(Sheet)
            xlwrite(FileOut,xlsout,Sheet);
        else
            xlwrite(FileOut,xlsout);
        end
    end
end
    
end