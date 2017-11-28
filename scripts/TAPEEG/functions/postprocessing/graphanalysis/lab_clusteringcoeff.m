% Calculate clustering coefficient
%
% cw = lab_clusteringcoeff(matrix)
%
%    Formula: Clustering coefficient = sum of all connections between
%             neighbours and connections to neighbours of N / sum of all
%             connections between neighbours of N
%
% Written by F. Hatz 2013

function cw = lab_clusteringcoeff(matrix)

nchans = size(matrix,1);
cw = zeros(nchans,1);

if size(matrix,1) ~= size(matrix,2)
    return
end

for i = 1:nchans
    ichans = setdiff(1:nchans,i);
    matrixtmp = tril(matrix(ichans,i) * matrix(i,ichans),-1);
    if sum(matrixtmp(:)) == 0
        cw(i,1) = 0;
    else
        matrixtmp2 = tril(matrix(ichans,ichans),-1) .* matrixtmp;
        cw(i,1) = sum(matrixtmp2(:)) / sum(matrixtmp(:));
    end
end

return