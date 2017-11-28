function matrixout = lab_rm_diagonal(matrixin)

if size(matrixin,1) == size(matrixin,2)
    chans = size(matrixin,1);
    diagonal = 1:chans+1:(chans*chans);
    include = setdiff((1:chans*chans),diagonal);
    matrixout = matrixin(include);
    matrixout = reshape(matrixout,chans-1,chans);
else
    disp('   error: removing diagonal not possible for nun square matrices')
end