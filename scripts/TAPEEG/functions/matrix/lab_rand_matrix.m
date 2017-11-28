function matrix = lab_rand_matrix(matrix)

[Nx,Ny] = size(matrix);
idx = find(tril(true(Nx,Ny),-1)==1);
Mvalues = matrix(idx);
Mdiag = diag(matrix);
matrix = zeros(Nx,Nx);
matrix(idx) = Mvalues(randperm(length(Mvalues)));
matrix = matrix + matrix';
matrix(1:Nx+1:end) = Mdiag;

end