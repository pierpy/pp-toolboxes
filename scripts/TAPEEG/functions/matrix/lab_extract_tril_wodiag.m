% Extract lower triangle of matrix without diagonal
%
% connections = lab_extract_tril_wodiag(matrix,labels)
%
% Written by F. Hatz 2013

function [connections,labels,X,Y] = lab_extract_tril_wodiag(matrix,labels)

if ~exist('labels','var')
    labels = {};
end
if isempty(matrix)
    connections = [];
    return
end

idx = tril(true(size(matrix(:,:,1))),-1);
connections = zeros(sum(idx(:)),size(matrix,3));
for i = 1:size(matrix,3)
    Mtmp = matrix(:,:,i);
    connections(:,i) = Mtmp(idx);
end

if ~isempty(labels) & nargout > 1
    labelstmp = {};
    for i = 1:length(labels)
        for j = i+1:length(labels)
            labelstmp{end+1,1} = [labels{i} '-' labels{j}]; %#ok<AGROW>
        end
    end
    labels = labelstmp;
end

if nargout > 2
    Y = mod(find(idx),size(matrix,1));
    Y(Y==0) = size(matrix,1);
    X = ceil(find(idx) / size(matrix,1));
else
    Y = [];
    X = [];
end

end