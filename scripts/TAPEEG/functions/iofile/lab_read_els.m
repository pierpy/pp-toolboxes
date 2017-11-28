% LOCS = lab_read_els(filename) - read Cartool *.els
% written by F. Hatz 2012

function LOCS = lab_read_els(filename)

% array = loadtxt(filename,'verbose','off');
array = lab_loadtxt(filename);
if strcmp(array{1,1},'ES01')
    LOCS.aux = 0;
    LOCS.auxlabels = {};
    array = array(4:end,:);
    m = 1;
    while m < size(array,1)
        if ischar(array{m,1})
            if array{m+2,1} == 0
                clustersize = 3 + array{m+1,1};
                LOCS.aux = LOCS.aux + array{m+1,1};
                if size(array,2) == 7 | size(array,2) == 8
                    LOCS.auxlabels = [LOCS.auxlabels array(m+3:m+2+array{m+1,1},7)'];
                else
                    LOCS.auxlabels = [LOCS.auxlabels array(m+3:m+2+array{m+1,1},4)'];
                end
                if m + clustersize < size(array,1)
                    
                    array = cat(1,array(1:m-1,:),array(m+clustersize:end,:));
                else
                    array = array(1:m-1,:);
                end
            else
                if m + 3 < size(array,1)
                    array = cat(1,array(1:m-1,:),array(m+3:end,:));
                else
                    array = array(1:m-1,:);
                end
            end
        end
        m = m + 1;
    end
    if size(array,2) == 7 | size(array,2) == 8
        LOCS.x = cell2mat(array(:,1))';
        LOCS.y = cell2mat(array(:,3))';
        LOCS.z = cell2mat(array(:,5))';
        LOCS.labels = array(:,7)';
    else
        LOCS.x = cell2mat(array(:,1))';
        LOCS.y = cell2mat(array(:,2))';
        LOCS.z = cell2mat(array(:,3))';
        LOCS.labels = array(:,4)';
    end
    for i = 1:size(LOCS.labels,2)
        if isnumeric(LOCS.labels{1,i})
            LOCS.labels{1,i} = num2str(LOCS.labels{1,i});
        else
            LOCS.labels{1,i} = regexprep(LOCS.labels{1,i},'_',' ');
        end
    end
    LOCS = lab_locs2sph(LOCS);
else
    LOCS = [];
    disp('    unkown els-format')
end