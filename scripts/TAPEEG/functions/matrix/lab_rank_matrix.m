% Rank matrix with selected exponent
%
% matrixout = lab_rank_matrix(matrix,exponent)
%
% written by F. Hatz

function matrixout = lab_rank_matrix(matrix,exponent)

if ~exist('exponent','var')
    exponent = 5;
end
disp(['      rank matrix with exponent ' num2str(exponent)])

for i = 1:size(matrix,3);
    matrixtmp = matrix(:,:,i);
    % store zero connections
    zmatrix = zeros(size(matrixtmp));
    zmatrix(matrixtmp == 0) = 1;
    
    % do ranking
    idx = find(tril(true(size(matrixtmp)),-1));
    tmp = matrixtmp(idx);
    [~,tmp2] = sort(tmp);
    tmp3 = (0:1/(length(tmp2)-1):1).^exponent;
    tmp3(tmp2) = tmp3;
    
    matrixtmp2 = zeros(size(matrixtmp));
    matrixtmp2(idx) = tmp3;
    matrixtmp2 = matrixtmp2 + matrixtmp2';
    matrixtmp2(zmatrix == 1) = 0;
    matrixout(:,:,i) = matrixtmp2;
end