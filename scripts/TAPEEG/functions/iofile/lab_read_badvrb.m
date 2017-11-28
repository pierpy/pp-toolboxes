function [badchans,filesize,bad] = lab_read_badvrb(file)

badchans = [];
bad.name = {};
bad.channels = {};

fid=fopen(file,'r');
if fid > 0
    tline=fgetl(fid);
    if strcmp(tline,'Detection of bad channels')
        tline=fgetl(fid);
        tline=fgetl(fid);
        tline=fgetl(fid);
        filesize = str2num(tline(13:end));
        tline=fgetl(fid);
        tline=fgetl(fid);
        while ~isnumeric(tline)
            bad.name{end+1} = tline(1:end-1);
            tline = fgetl(fid);
            bad.channels{end+1} = str2num(tline);
            badchans = union(badchans,bad.channels{end});
            tline=fgetl(fid);
            tline=fgetl(fid);
        end
        fclose(fid);
    else
        badchans = -1;
        filesize = 0;
    end
else
    badchans = -1;
    filesize = 0;
end
clearvars tline fid;

if size(badchans,1) > 1
    badchans = badchans';
end

return