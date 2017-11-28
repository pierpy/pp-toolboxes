% Calculate Minimum Spanning Tree
%
% [MST,MSTweight,MSTline] = lab_MST(matrix)
%
% written by F.Hatz 2014

function [MST,MSTweight,MSTline] = lab_MST(matrix)

% calculate inverse matrix
matrixI = abs(matrix).^-1;
matrixI(1:(size(matrix,1)+1):end) = 0;

tmp = diag(matrixI);
matrixI = tril(matrixI) + tril(matrixI)';
matrixI(1:size(matrixI,1)+1:end) = tmp;
clearvars tmp

MST = full(kruskal_mst(sparse(matrixI)));
MSTweight = sum(MST(:));
MST(MST>0) = 1;

[tmp1,tmp2,tmp3] = find(MST);
MSTline = [tmp1 tmp2 tmp3];
clearvars tmp1 tmp2 tmp3