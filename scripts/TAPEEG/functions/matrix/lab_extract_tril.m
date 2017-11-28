% Extract lower triangle of matrix
%
% connections = lab_extract_tril(matrix)
%
% Written by F. Hatz 2013

function [connections,labels] = lab_extract_tril(matrix,labels,includediag)

if ~exist('includediag','var')
    includediag = true;
end
if ~exist('labels','var')
    labels = {};
end
if isempty(matrix)
    connections = [];
    return
end

if max(max(abs(matrix(:,:,1) - matrix(:,:,1)'))) > 10^-10
    issymm = false;
else
    issymm = true;
end

if issymm == false
    idx = true(size(matrix,1),size(matrix,2));
elseif includediag == true
    idx = tril(true(size(matrix(:,:,1))));
else
    idx = tril(true(size(matrix(:,:,1))),-1);
end
connections = zeros(sum(idx(:)),size(matrix,3));
for i = 1:size(matrix,3)
    Mtmp = matrix(:,:,i);
    connections(:,i) = Mtmp(idx);
end

if ~isempty(labels)
    labelstmp = {};
    if issymm == false
        for i = 1:length(labels)
            for j = 1:length(labels)
                if includediag == true | i ~= j
                    labelstmp{end+1,1} = [labels{i} '-' labels{j}]; %#ok<AGROW>
                end
            end
        end
    else
        for i = 1:length(labels)
            for j = i:length(labels)
                if includediag == true | i ~= j
                    labelstmp{end+1,1} = [labels{i} '-' labels{j}]; %#ok<AGROW>
                end
            end
        end
    end
    labels = labelstmp;
end