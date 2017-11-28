% Helper function for lab_collect_spectras
%
% [calc,result] = lab_collect_spectraldata_relativepower(cfg,calc,result,SpectAllMean,SpectAllMedian)
%
% written by F. Hatz 2012

function [calc,result] = lab_collect_spectraldata_relativepower(cfg,calc,result,SpectAllMean,SpectAllMedian)

% Collect realtive power of every channel
if isfield(cfg.CollectFFT,'source') & strcmp(cfg.CollectFFT.source,'median')
    if ~isempty(SpectAllMedian)
        SpectAllM2 = SpectAllMedian;
    else
        SpectAllM2 = [];
    end
else
    if ~isempty(SpectAllMean)
        SpectAllM2 = SpectAllMean;
    else
        SpectAllM2 = [];
    end
end
if ~isempty(SpectAllM2)
    for j = 1:size(SpectAllM2,1);
        SpectAllM2(j,:) = SpectAllM2(j,:) ./ sum(SpectAllM2(j,:));
    end
    if strcmp('empty',result.FBM)
        result.FBM = SpectAllM2;
        result.FBMpat = cellstr(calc.patient);
    else
        result.FBM = cat(3,result.FBM,SpectAllM2);
        result.FBMpat = [result.FBMpat cellstr(calc.patient)];
    end
end
clearvars SpectAllM2