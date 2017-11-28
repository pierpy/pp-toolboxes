function [Filename,Filepath,Format,FilenameS,Segment] = lab_filename(Filename)

tmp = strfind(Filename,filesep);
if ~isempty(tmp) & tmp(end) == length(Filename)
    if length(tmp) == 1
        tmp = [];
    else
        tmp = tmp(1:end-1);
    end
end
if ~isempty(tmp)
    Filepath = Filename(1:tmp(end));
    Filename = (Filename(tmp(end)+1:end));
else
    Filepath = [];
end

tmp = strfind(Filename,'.');
if ~isempty(tmp) & tmp(end) >= (length(Filename)-5)
    Format = lower(Filename(tmp(end)+1:end));
    FilenameS = Filename(1:tmp(end)-1);
else
    FilenameS = Filename;
    Format = '';
end

tmp = strfind(FilenameS,'_');
if ~isempty(tmp) & tmp(end) == length(FilenameS)
    if length(tmp) > 1
        tmp = tmp(1:end-1);
    else
        tmp = [];
    end
end
for i = 1:size(tmp,2)
    if strcmp(FilenameS(tmp(i)+1),'S')
        Sstart = tmp(i) + 2;
        Ssize = 0;
        while (Sstart+Ssize) <= length(FilenameS) & ~isempty(str2num(FilenameS(Sstart + Ssize))) %#ok<ST2NM>
            Ssize = Ssize + 1;
        end
        if Ssize > 0
            Segment = str2num(FilenameS(Sstart:(Sstart+Ssize-1))); %#ok<ST2NM>
        end
    end
end
if ~exist('Segment','var')
    Segment = [];
end
if ~isempty(FilenameS)
    FilenameS = regexprep(FilenameS,':','_');
    FilenameS = regexprep(FilenameS,';','_');
end