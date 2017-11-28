% Calculate average clustering coefficient and average apth length for random
% matrices
%
% [cw,lambda] = lab_graph_randmatrix(matrix,Rnumber,Riter,docw,dolambda)
%
% Written by F. Hatz 2013

function [cw,lambda] = lab_graph_randmatrix(matrix,Rnumber,Riter,docw,dolambda)

if~exist('docw','var')
    docw = 1;
end
if~exist('dolambda','var')
    dolambda = 1;
end

rng('default');
rng('shuffle');
[Nx,Ny] = size(matrix);

Ndots = floor(Rnumber/10);
if Nx == Ny
    fprintf(['RandMatrix(' num2str(Rnumber) ') '])
    lambda = zeros(1,Rnumber);
    cw = zeros(1,Rnumber);
    for Rnum = 1:Rnumber
        matrixtmp = lab_rand_matrix_fixed(matrix,Riter);
        matrixI = abs(matrixtmp).^-1;
        matrixI(1:(size(matrixtmp,1)+1):end) = 0;
        
        if dolambda == 1
            % Path length
            distance = (distance_wei(matrixI)).^-1;
            distance(1:Nx+1:end) = 0;
            lambda(Rnum) = (sum(distance(:)) / (Nx*(Nx-1))).^-1;
        end
        
        if docw == 1
            % Clustering coeff, Mika Rubinov, UNSW, 2007-2010
            K=sum(matrixtmp~=0,2);
            cyc3=diag((matrixtmp.^(1/3))^3);
            K(cyc3==0)=inf;             %if no 3-cycles exist, make C=0 (via K=inf)
            cw(Rnum)=mean(cyc3./(K.*(K-1)));
        end
        if mod(Rnum,Ndots) == 0
            fprintf('.')
        end
    end
    cw = mean(cw);
    lambda = mean(lambda);
else
    disp('    error: input matrix must be squared')
end

end