% Helper function for lab_collect_spectras
%
% [cogfreq,Tcogfreq] = lab_collect_spectraldata_mf(SpectAllM,SpectAllF,settings)
%
% written by F. Hatz 2012

function [cogfreq,Tcogfreq] = lab_collect_spectraldata_mf(SpectAllM,SpectAllF,settings)

if ~isnumeric(SpectAllF)
    cogfreq = NaN;
    Tcogfreq = [0 0];
    return
end

freqlow = find(SpectAllF >= settings.lowfreqcog, 1 );
freqhigh = find(SpectAllF <= settings.highfreqcog, 1,'last');
if ~isempty(freqlow) & ~isempty(freqhigh)
    Spect2 = SpectAllM(:,freqlow:freqhigh);
    [rc,cc] = ndgrid(1:size(Spect2,1),1:size(Spect2,2));
    Spectt = sum(Spect2(:));
    c1 = sum(Spect2(:) .* rc(:)) / Spectt;
    c2 = sum(Spect2(:) .* cc(:)) / Spectt;
    cogfreq = SpectAllF(freqlow -1) + (c2 * (SpectAllF(2) - SpectAllF(1)));
    tmp = ((cogfreq - SpectAllF(1))/(SpectAllF(2)-SpectAllF(1)))+1;
    Tcogfreq = [tmp tmp];
else
    cogfreq = NaN;
    Tcogfreq = [0 0];
end