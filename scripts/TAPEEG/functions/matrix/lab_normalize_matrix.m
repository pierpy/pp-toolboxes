% Normalize matrix (minimal value = 0 / maximal value = 1)
%
% matrix = lab_normalize_matrix(matrix,keepdiag)
%      keepdiagonal: true = diagonal remains unchanged
%
% written by F. Hatz

function matrix = lab_normalize_matrix(matrix,keepdiag)

if ~exist('keepdiag','var')
    keepdiag = false;
end

if keepdiag == false
    numchans = size(matrix,1);
    for i = 1:size(matrix,3);
        tmp = matrix(:,:,i);
        tmp(1:numchans+1:end) = NaN;
        matrix(:,:,i) = tmp;
    end
end

matrix = (matrix - min(matrix(:))) / (max(matrix(:)) - min(matrix(:)));

if keepdiag == false
    for i = 1:size(matrix,3);
        tmp = matrix(:,:,i);
        tmp(1:numchans+1:end) = 0;
        matrix(:,:,i) = tmp;
    end
end
 