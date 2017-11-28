% Read electrodes localications in .xyz format
%
% LOCS = lab_read_xyz(filename)
%
% written by F. Hatz 2012

function LOCS = lab_read_xyz(filename)

% array = loadtxt(filename,'verbose','off');
array = lab_loadtxt(filename);
if size(array,2) == 6
    LOCS.x = cell2mat(array(2:end,1))';
    LOCS.y = cell2mat(array(2:end,3))';
    LOCS.z = cell2mat(array(2:end,5))';
    LOCS.labels = array(2:end,6)';
elseif size(array,2) == 7
    LOCS.x = cell2mat(array(2:end,1))';
    LOCS.y = cell2mat(array(2:end,3))';
    LOCS.z = cell2mat(array(2:end,5))';
    LOCS.labels = array(2:end,7)';
elseif isempty(array{1,3})
    numelec = array{1,1};
    LOCS.x = cell2mat(array(2:numelec+1,1))';
    LOCS.y = cell2mat(array(2:numelec+1,2))';
    LOCS.z = cell2mat(array(2:numelec+1,3))';
    LOCS.labels = array(2:numelec+1,4)';
else
    LOCS.x = cell2mat(array(:,1))';
    LOCS.y = cell2mat(array(:,2))';
    LOCS.z = cell2mat(array(:,3))';
    LOCS.labels = array(:,4)';
end
for i = 1:size(LOCS.labels,2)
    if isnumeric(LOCS.labels{1,i})
        LOCS.labels{1,i} = num2str(LOCS.labels{1,i});
    end
end

LOCS = lab_locs2sph(LOCS);

return