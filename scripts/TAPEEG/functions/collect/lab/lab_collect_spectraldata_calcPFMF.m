% Helper function for lab_collect_spectras
%
% [R,T] = lab_collect_spectraldata_calcPFMF(SpectAllM,SpectAllF,settings,header)
%
% written by F. Hatz 2012

function [R,T] = lab_collect_spectraldata_calcPFMF(SpectAllM,SpectAllF,settings,header)

if ~exist('header','var')
    header = [];
end

% Calculate background peak / median freq
if settings.calcsingle == false
    if settings.domedian == true
        SpectAllM2 = median(SpectAllM(settings.mappingBA,:),1);
    else
        SpectAllM2 = mean(SpectAllM(settings.mappingBA,:).^0.5,1).^2;
    end
else
    SpectAllM2 = SpectAllM(settings.mappingBA,:);
end
if ~isfield(settings,'PF') & isnumeric(SpectAllF)
    R = lab_collect_spectraldata_quality(SpectAllM2,SpectAllF,settings,header);
    if ~isnan(R.peakfreq)
        settings.PF = find(SpectAllF >= R.peakfreq, 1 );
    end
end
[R,T] = lab_collect_spectraldata_quality(SpectAllM2,SpectAllF,settings,header);
T.mappingBA = settings.mappingBA;

% Calculate for every channel
Rtmp = lab_collect_spectraldata_quality(SpectAllM,SpectAllF,settings,header);
R.peakfreqAll = Rtmp.peakfreqAll;
R.cogfreqAll = Rtmp.cogfreqAll;
clearvars Rtmp

% Calculate Mappings
if isfield(settings,'mappings') & ~isempty(settings.mappings)
    for j = 1:size(settings.mappings,2)
        if settings.calcsingle == false
            if settings.domedian == false
                SpectAllM2 = mean(SpectAllM(settings.mappings{1,j},:).^0.5,1).^2;
            else
                SpectAllM2 = median(SpectAllM(settings.mappings{1,j},:),1);
            end
            Rtmp = lab_collect_spectraldata_quality(SpectAllM2,SpectAllF,settings,header);
            R.cogfreqMap(j,1) = Rtmp.cogfreq;
            R.peakfreqMap(j,1) = Rtmp.peakfreq;
            clearvars Rtmp SpectAllM2
        elseif settings.domedian == true
            R.cogfreqMap(j,1) = lab_collect_spectraldata_mf(SpectAllM(settings.mappings{1,j},:),SpectAllF,settings);
            R.peakfreqMap(j,1) = median(R.peakfreqAll(settings.mappings{1,j}));
        end
    end
    T.mappings = settings.mappings;
end

T.SpectAll = SpectAllM;
if isfield(settings,'Valid')
    T.Valid = settings.Valid;
end