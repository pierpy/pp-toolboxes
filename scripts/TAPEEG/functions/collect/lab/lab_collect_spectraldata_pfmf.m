% Helper function for lab_collect_spectras
%
% [cfg,result] = lab_collect_spectraldata_pfmf(cfg,result,calc,SpectAllF,SpectAllMean,SpectAllMedian,Mappings,Valid,header)
%
% written by F. Hatz 2012

function [cfg,result] = lab_collect_spectraldata_pfmf(cfg,result,calc,SpectAllF,SpectAllMean,SpectAllMedian,Mappings,Valid,header)

settings = cfg.CollectFFT;
settings.mappingBA = calc.mappingBA;
if isfield(cfg.CollectFFT,'source') & strcmp(cfg.CollectFFT.source,'median') & ~isempty(SpectAllMedian)
    settings.domedian = true;
else
    settings.domedian = false;
    cfg.CollectFFT.source = 'mean';
end
if ~isfield(settings,'calcsingle')
    settings.calcsingle = false;
end
if exist('Valid','var')
    settings.Valid = Valid;
else
    settings.Valid = [];
end
if exist('Mappings','var') & isfield(Mappings,'mappings')
    settings.mappings = Mappings.mappings;
else
    settings.mappings = [];
end

if settings.domedian == true
    if ~isempty(SpectAllMedian)
        SpectAllM = SpectAllMedian;
    else
        SpectAllM = [];
    end
else
    if ~isempty(SpectAllMean)
        SpectAllM = SpectAllMean;
    else
        SpectAllM = [];
    end
end
if isempty(SpectAllM)
    disp(['    Calculation of peak & median freq for ' calc.patient ' not possible'])
    return
end

[R,T] = lab_collect_spectraldata_calcPFMF(SpectAllM,SpectAllF,settings,header);

if ~isfield(result,'T')
    result.T = T;
    result.R = R;
    result.patient = cellstr(calc.patient);
else
    result.T = [result.T T];
    result.R = [result.R R];
    result.patient = [result.patient cellstr(calc.patient)];
end

end