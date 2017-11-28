function [ meanDur, meanOcc, gev, totalGev, meanCov, statesSequence, singleCorrelations ] = extractbackfittingstruct( self, backfittingResults )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
    meanDur = backfittingResults.mean_dur;
    meanOcc = backfittingResults.freq;
    gev = backfittingResults.gev;
    totalGev = backfittingResults.GEV_ALL;
    meanCov = backfittingResults.cov;
    statesSequence = backfittingResults.states_sequence;
    singleCorrelations = backfittingResults.single_epoch_corr;
end

