% Extract upper triangle of matrix without diagonal
%
% connections = lab_extract_triu_wodiag(matrix,labels)
%
% Written by F. Hatz 2013

function [connections,labels,X,Y] = lab_extract_triu_wodiag(matrix,labels)

if ~exist('labels','var')
    labels = {};
end
if isempty(matrix)
    connections = [];
    return
end

idx = triu(true(size(matrix(:,:,1))),1);
connections = zeros(sum(idx(:)),size(matrix,3));
for i = 1:size(matrix,3)
    Mtmp = matrix(:,:,i);
    connections(:,i) = Mtmp(idx);
end

if ~isempty(labels) & nargout > 1
    labelstmp = {};
    for i = 2:length(labels)
        for j = 1:i-1
            labelstmp{end+1,1} = [labels{i} '-' labels{j}]; %#ok<AGROW>
        end
    end
    labels = labelstmp;
end

if nargout > 2
    X = mod(find(idx),size(matrix,1));
    X(X==0) = size(matrix,1);
    Y = ceil(find(idx) / size(matrix,1));
else
    Y = [];
    X = [];
end

end