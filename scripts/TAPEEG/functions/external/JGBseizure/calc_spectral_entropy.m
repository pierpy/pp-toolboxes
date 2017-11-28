function S = calc_spectral_entropy(EEG_power,f)

% referentie: {Viertiö-Oja, 2004 #33}
P = EEG_power(find(f>0.5,1):end);
Pn = P/sum(P);

S = sum(Pn.*log(1./Pn))/log(length(Pn));




