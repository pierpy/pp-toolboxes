function [matrix,eff] = lab_rand_matrix_fixed(matrix,ITER)

n=size(matrix,1);
Diag = diag(matrix);
matrix(1:n+1:end) = 0;
idx = tril(true(size(matrix)),-1);
[i,j]=find(idx);
K = sum(idx(:));
ITER=K*ITER;

matrix = tril(matrix,-1) + tril(matrix,-1)';
limit = (max(matrix(:)) - min(matrix(:))) / 10;
eff = 0;
maxAttempts= round(n*K/(n*(n-1)));
for N = 1:ITER
    att=0;
    while (att<=maxAttempts) %while not rewired
        while 1
            e1=ceil(K*rand);
            e2=ceil(K*rand);
            while (e2==e1),
                e2=ceil(K*rand);
            end
            a=i(e1); b=j(e1);
            c=i(e2); d=j(e2);

            if all(a~=[c d]) && all(b~=[c d]);
                break %all four vertices must be different
            end
        end

        if rand>0.5
            i(e2)=d; j(e2)=c; 	%flip edge c-d with 50% probability
            c=i(e2); d=j(e2); 	%to explore all potential rewirings
        end
        
        %rewiring condition
        if abs(matrix(a,b)-matrix(c,d)) < limit & abs(matrix(a,c)-matrix(b,d)) < limit
            matrix(a,b) = matrix(b,d);
            matrix(b,d) = matrix(d,c);
            matrix(d,c) = matrix(c,a);
            matrix(c,a) = matrix(b,a);
            matrix(b,a) = matrix(a,b);
            matrix(d,b) = matrix(b,d);
            matrix(c,d) = matrix(d,c);
            matrix(a,c) = matrix(c,a);
            
            j(e1) = d; %reassign edge indices
            j(e2) = b;
            eff = eff+1;
            break;
        end %rewiring condition
        att=att+1;
    end
end
matrix(1:n+1:end) = Diag;
