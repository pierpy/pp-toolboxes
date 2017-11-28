% Read electrodes localications in .sfp format
%
% LOCS = lab_read_sfp(filename)
%
% written by F. Hatz 2012

function LOCS = lab_read_sfp(filename)

% array = loadtxt(filename,'verbose','off');
array = lab_loadtxt(filename);
LOCS.x = cell2mat(array(:,3))';
LOCS.y = -cell2mat(array(:,2))';
LOCS.z = cell2mat(array(:,4))';
LOCS.labels = array(:,1)';
for i = 1:size(LOCS.labels,2)
    if isnumeric(LOCS.labels{1,i})
        LOCS.labels{1,i} = num2str(LOCS.labels{1,i});
    end
end

LOCS = lab_locs2sph(LOCS);

return