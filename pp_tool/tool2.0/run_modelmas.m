function run_modelmas(confstruct, eeg, seed)
s = RandStream('mt19937ar','Seed',seed);
RandStream.setGlobalStream(s);

[eeg, gfp_peak_indices, gfp_curve] = modmaps_preprocess(confstruct, eeg);