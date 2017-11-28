% Read electrodes localications in .sxyz format
%
% LOCS = lab_read_sxyz(filename)
%
% written by F. Hatz 2012

function LOCS = lab_read_sxyz(filename)

% array = loadtxt(filename,'verbose','off');
array = lab_loadtxt(filename);
LOCS.x = cell2mat(array(2:end,1))';
LOCS.y = cell2mat(array(2:end,2))';
LOCS.z = cell2mat(array(2:end,3))';
LOCS.labels = array(2:end,4)';
for i = 1:size(LOCS.labels,2)
    if isnumeric(LOCS.labels{1,i})
        LOCS.labels{1,i} = num2str(LOCS.labels{1,i});
    end
end
LOCS = lab_locs2sph(LOCS);

return