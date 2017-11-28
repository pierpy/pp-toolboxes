function [Filefirst,Filelast] = lab_find_lastperfolder(Filelist)

for i = 1:length(Filelist);
    [~,Filepath{i}] = lab_filename(Filelist{i}); %#ok<AGROW>
end
[~,tmp] = unique(Filepath,'last');
Filelast = zeros(1,length(Filelist));
Filelast(tmp) = 1;
[~,tmp] = unique(Filepath,'first');
Filefirst = zeros(1,length(Filelist));
Filefirst(tmp) = 1;

end