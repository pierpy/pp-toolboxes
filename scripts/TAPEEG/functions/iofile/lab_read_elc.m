% LOCS = lab_read_elc(filename) - read *.elc
% written by F. Hatz 2014

function LOCS = lab_read_elc(filename)

% array = loadtxt(filename,'verbose','off');
array = lab_loadtxt(filename);

if ~strcmp(array{1,1},'#') | ~strcmp(array{1,2},'ASA')
    return
else
    array = array(2:end,:);
end
if strcmp(array{1,1},'ReferenceLabel')
    LOCS.reference = array{1,2};
    array = array(2:end,:);
end
if strcmp(array{1,1},'UnitPosition')
    if strcmp(array{1,2},'mm')
        Factor = 10;
    elseif strcmp(array{1,2},'m')
        Factor = 0.1;
    else
        Factor = 1;
    end
    array = array(2:end,:);
else
    Factor = 1;
end
if strcmp(array{1,1},'NumberPositions=')
    Nchannels = array{1,2};
    array = array(2:end,:);
else
    Nchannels = [];
end
while ischar(array{1,1})
    if size(array,1) > 1
        array = array(2:end,:);
    else
        disp('    unkown elc-format')
        LOCS = [];
        return
    end
end
if isempty(Nchannels) & 2*Nchannels ~= size(array,1)-1
    Nchannels = find(cellfun(@isnumeric,array(:,1)),1,'last');
end
if strcmp(array{Nchannels+1,1},'Labels') | size(array,1) == 2*Nchannels+1
    LOCS.labels = array(Nchannels+2:end,1)';
else
    LOCS.labels = cellstr(num2str((1:Nchannels)'))';
end
try
    LOCS.x = cell2mat(array(1:Nchannels,1))' / Factor;
    LOCS.y = cell2mat(array(1:Nchannels,2))' / Factor;
    LOCS.z = cell2mat(array(1:Nchannels,3))' / Factor;
    LOCS = lab_locs2sph(LOCS);
catch %#ok<CTCH>
    disp('    unkown elc-format')
    LOCS = [];
end

end