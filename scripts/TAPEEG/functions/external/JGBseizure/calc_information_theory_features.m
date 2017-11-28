function information_theory_features = calc_information_theory_features(epoch,EEG_power,f)

%Shannon entropy
Epoch_entropy = entropy(epoch');

%singular value decomposition entropy en Fisher information
SVD_entropy = calc_svd_entropy_and_fisher_info(epoch);

%spectral Entropy
 Spectral_entropy = calc_spectral_entropy(EEG_power,f);

%collect information theory features
 information_theory_features = [Epoch_entropy, SVD_entropy,...
     fisher_information, Spectral_entropy];
