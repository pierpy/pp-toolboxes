function result = lab_get_group(subjects,results)

if ~exist('results','var')
    [tmp1,tmp2] = uigetfile('*.xlsx','Select file with results');
    results = fullfile(tmp2,tmp1);
    results = lab_read_xls(results);
end

for i = 1:length(subjects)
    tmp = strfind(subjects{i},'_');
    if ~isempty(tmp)
        subjects{i} = subjects{i}(1:tmp(1)-1);
    end
end

if isempty(results)
    subjectsR = {''};
elseif  min(cellfun(@ischar,results(1,:))) == 1
    tmp = intersect(results(1,:),subjects);
    if ~isempty(tmp) & size(results,1) > 1
        subjectsR = results(1,:)';
        results = results(2,:)';
    end
elseif  size(results,2) > 1 & min(cellfun(@ischar,results(1,2:end))) == 1
    tmp = intersect(results(1,2:end),subjects);
    if ~isempty(tmp) & size(results,1) > 1
        subjectsR = results(1,2:end)';
        results = results(2,2:end)';
    end
elseif  min(cellfun(@ischar,results(:,1))) == 1
    tmp = intersect(results(:,1),subjects);
    if ~isempty(tmp) & size(results,2) > 1
        subjectsR = results(:,1)';
        results = results(:,2)';
    end
elseif  size(results,2) > 1 & min(cellfun(@ischar,results(2:end,1))) == 1
    tmp = intersect(results(2:end,1),subjects);
    if ~isempty(tmp) & size(results,2) > 1
        subjectsR = results(2:end,1)';
        results = results(2:end,2)';
    end
else
    subjectsR = {''};
end

for i = 1:length(subjects)
    tmp = find(strcmp(subjectsR,subjects{i}));
    if ~isempty(tmp)
        result(1,i) = results{tmp,1};
    else
        result(1,i) = 0;
    end
end