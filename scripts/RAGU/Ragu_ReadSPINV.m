function [res] = Ragu_ReadSPINV(infile)

fid = fopen(infile,'rb');
if fid == -1
    res = [];
    return;
end

f = fread(fid,'single');
fclose(fid);

res.nvox = f(1);
res.nchan = f(2);
res.spinv = reshape(f(3:(res.nvox*res.nchan*3+2)),[res.nchan,3,res.nvox]);
