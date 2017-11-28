function matrixout = lab_add_diagonal(matrixin)

if (size(matrixin,1)+1) == size(matrixin,2)
    chans = size(matrixin,2);
    diagonal = 1:chans+1:(chans*chans);
    include = setdiff((1:chans*chans),diagonal);
    matrixout = zeros(chans,chans);
    matrixout(include) = matrixin;
else
    disp('   error: adding diagonal only possible for matrixes(x-1*x)')
end