function [ expVar, prototypes, statesSequence, cv, kl ] = extractsegmentationstruct( self, clusteringResults, ntemplatesToExtract )
% extractsegmentationstruct: returns segementations results each in a
% single variables.
%   Detailed explanation goes here
    expVar = clusteringResults.exp_var;
    prototypes = clusteringResults.template{ntemplatesToExtract};
    statesSequence = clusteringResults.clusters{ntemplatesToExtract};
    cv = clusteringResults.CV;
    kl = clusteringResults.kl;
end

