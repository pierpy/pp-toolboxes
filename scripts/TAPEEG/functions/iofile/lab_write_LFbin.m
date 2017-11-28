% Write Cartool leadfield .bin
%
% lab_write_LFbin(filename,LFbin)
%
% written by F. Hatz 2012

function lab_write_LFbin(filename,LFbin)

LFbin = permute(LFbin,[2 3 1]);

fid=fopen(filename,'w');
fwrite(fid,LFbin,'float32');
fclose(fid);