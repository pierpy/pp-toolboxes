function pos = ReadXYZ(fname)

% Copyright 2009-2011 Thomas Koenig
% distributed under the terms of the GNU AFFERO General Public License

[fid,message] = fopen(fname,'rt');

if fid == -1
    errordlg(message,'Read XYZ File');
    pos = [];
    return;
end

l = fgets(fid);
[params,pcount] = sscanf(l,'%f')

C = textscan(fid,'%f%f%f%s');

fclose(fid);

x = C{1};
y = C{2};
z = C{3};

pos = [-x y z];

if pcount > 1
    pos = pos / params(2);
else
    nrm = max(sqrt(sum(pos.^2,2)));
    pos = pos / nrm;
end
%pos = [y -x z] / params(2);
