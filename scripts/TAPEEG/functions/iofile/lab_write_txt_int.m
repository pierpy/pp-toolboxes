% Write eeg/meg in ascii (txt)
%
% lab_write_txt(filename,data,header)
%
% written by F. Hatz 2012

function lab_write_txt_int(filename,data,header)
[~,filepath,~,filename] = lab_filename(filename);
filename = fullfile(filepath,[filename '.txt']);

fid=fopen(filename,'w');
for i = 1:size(data,1)
    fprintf(fid,header.channels(i,:));
    fprintf(fid,'\t');
end
fprintf(fid,native2unicode([13 10]));
fclose(fid);
dlmwrite(filename,data','delimiter','\t','newline','pc','-append');

return
