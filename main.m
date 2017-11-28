clear
confstruct.tf = 15000;
confstruct.fs = 250;
confstruct.nch = 19;
confstruct.minclusters = 3;
confstruct.maxclusters = 6;
data = randn(15000, 19);
MM = mustates (confstruct);

[ eeg, gfp_peak_indices, gfp_curve ] = MM.preprocess( data );
[ clusteringResults ] = MM.segmentation( data, gfp_peak_indices );
[ expVar, prototypes, statesSequence, cv, kl ] = MM.extractsegmentationstruct(clusteringResults, 4 ); 
[ backfittingResults ] = MM.computemsparameters( eeg, prototypes);
