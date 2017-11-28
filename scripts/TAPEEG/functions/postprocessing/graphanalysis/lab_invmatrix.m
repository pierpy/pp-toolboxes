% Inverse matrix
%
% [matrixI,matrixlineI] = lab_invmatrix(matrix)
%
% Written by F. Hatz 2013

function [matrixI,matrixlineI] = lab_invmatrix(matrix)

matrixI = abs(matrix).^-1;
matrixI(1:(size(matrix,1)+1):end) = 0;
[tmp1,tmp2,tmp3] = find(matrixI);
matrixlineI = [tmp1 tmp2 tmp3];
clearvars tmp1 tmp2 tmp3

% Old Formula, not used anymore
%
% if min(matrix(:)) >= 0 & max(matrix(:)) <= 1
%     matrixI = -abs(matrix);
%     matrixI(matrixI == 0) = NaN;
%     matrixI(1:(size(matrix,1)+1):end) = NaN;
%     matrixI = (matrixI - min(matrixI(:)));
%     matrixI = matrixI / max(matrixI(:));
%     matrixI(1:(size(matrix,1)+1):end) = 0;
%     matrixI(isnan(matrixI)) = Inf;
%     matrixI(1:size(matrixI,1)+1:end) = 0;
% else
%     matrixI = abs(matrix).^-1;
%     matrixI(1:(size(matrix,1)+1):end) = 0;
% end

