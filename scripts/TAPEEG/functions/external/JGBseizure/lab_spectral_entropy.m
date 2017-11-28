function S = lab_spectral_entropy(P)

Pn = P/sum(P);
S = sum(Pn.*log(1./Pn))/log(length(Pn));

end