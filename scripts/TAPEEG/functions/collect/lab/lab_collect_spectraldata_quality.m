% Helper function for lab_collect_spectras
%
% [R,T] = lab_collect_spectraldata_quality(SpectAllM,SpectAllF,settings,header)
%
% written by F. Hatz 2012

function [R,T] = lab_collect_spectraldata_quality(SpectAllM,SpectAllF,settings,header)

if ~exist('header','var')
    header = [];
end
if ~exist('settings','var')
    settings = [];
end
if ~isfield(settings,'domedian')
    settings.domedian = [];
end
if ~isfield(settings,'PF')
    settings.PF = [];
end
if nargout == 2
    doquality = true;
else
    doquality = false;
end

numchans = size(SpectAllM,1);
if isempty(numchans)
    numchans = 1;
end

R.cogfreq = NaN(numchans,1);
R.cogfreqAll = NaN(numchans,1);
R.cogareapower = zeros(numchans,1);
R.peakfreq = NaN(numchans,1);
R.peakfreqAll = NaN(numchans,1);
R.peakamp = zeros(numchans,1);
R.peakratio = zeros(numchans,1);
R.areapower = zeros(numchans,1);
R.peak2min = zeros(numchans,1);
R.bplast = zeros(numchans,1);
T.cogfreq = zeros(numchans,2);
T.cogareapower = zeros(numchans,2);
T.peakfreq = zeros(numchans,2);
T.peakamp = zeros(numchans,2);
T.peakratio = zeros(numchans,2);
T.areapower = zeros(numchans,2);
T.peak2min = zeros(numchans,2);
T.bplast = zeros(numchans,2);
T.SpectF = SpectAllF;
T.domedian = settings.domedian;

if isempty(SpectAllM)
    return
end
if ~isnumeric(SpectAllF)
    if settings.domedian == true
        T.Spect = median(SpectAllM,1);
    else
        T.Spect = mean(SpectAllM .^0.5,1).^2;
    end
    return
end

if ~isfield(settings,'MinPeak2Min') | isempty(settings.MinPeak2Min)
    settings.MinPeak2Min = 1.3;
end
if ~isfield(settings,'lowfreqpeak')
    settings.lowfreqpeak = SpectAllF(1);
    settings.highfreqpeak = SpectAllF(end);
end
if ~isfield(settings,'lowfreqcog')
    settings.lowfreqcog = SpectAllF(2);
    settings.highfreqcog = SpectAllF(end);
end

for Nspect = 1:numchans
    Spect = SpectAllM(Nspect,:);
    if isfield(settings,'PF') & ~isempty(settings.PF) & ~isnan(settings.PF)
        freqlow = round(settings.PF - 2*(settings.qualityrange / (SpectAllF(2)-SpectAllF(1))));
        if ~isempty(settings.lowfreqpeak) & freqlow < find(SpectAllF >= settings.lowfreqpeak,1) & ...
                settings.PF >= find(SpectAllF >= settings.lowfreqpeak,1)
            freqlow = find(SpectAllF >= settings.lowfreqpeak, 1 );
        elseif freqlow < 1
            freqlow = 1;
        end
        freqhigh = round(settings.PF + 2*(settings.qualityrange / (SpectAllF(2)-SpectAllF(1))));
        if isfield(settings,'highfreqpeak') & ~isempty(settings.highfreqpeak) & freqhigh > find(SpectAllF <= settings.highfreqpeak, 1,'last') & ...
                settings.PF <= find(SpectAllF <= settings.highfreqpeak, 1,'last')
            freqhigh = find(SpectAllF <= settings.highfreqpeak, 1,'last');
        elseif freqhigh > length(SpectAllF)
            freqhigh = length(SpectAllF);
        end
    elseif isfield(settings,'PF') & ~isempty(settings.PF) & isnan(settings.PF)
        freqlow = [];
        freqhigh = [];
    else
        freqlow = find(SpectAllF >= settings.lowfreqpeak,1);
        freqhigh = find(SpectAllF <= settings.highfreqpeak,1,'last');
    end
    if ~isempty(freqlow) & ~isempty(freqhigh)
        [R.peakfreqAll(Nspect,1),T.peakfreq(Nspect,:)] = calcPF(Spect,SpectAllF,freqlow,freqhigh,settings.MinPeak2Min);
    end
    freqlow = find(SpectAllF >= settings.lowfreqpeak, 1 );
    freqhigh = find(SpectAllF <= settings.highfreqpeak, 1,'last');
    Spect2 = Spect(1,freqlow:freqhigh);
    R.areapower(Nspect,1) = trapz(SpectAllF(1,freqlow:freqhigh),Spect2,2);
    T.areapower(Nspect,:) = [freqlow freqhigh];
    if ~isnan(R.peakfreqAll(Nspect,1)) & doquality == true
        % calculate quality check
        peaklow = find(SpectAllF >= R.peakfreqAll(Nspect,1) - settings.qualityrange, 1 );
        peakhigh = find(SpectAllF <= R.peakfreqAll(Nspect,1) + settings.qualityrange, 1,'last');
        if peaklow < peakhigh
            R.peakamp(Nspect,1) = trapz(SpectAllF(1,peaklow:peakhigh),Spect(1,peaklow:peakhigh),2);
            T.peakamp(Nspect,:) = [peaklow peakhigh];
            R.peakratio(Nspect,1) = trapz(SpectAllF(1,peaklow:peakhigh),Spect(1,peaklow:peakhigh),2) / ...
                trapz(SpectAllF(1,freqlow:freqhigh),Spect2,2);
            T.peakratio(Nspect,:) = [freqlow freqhigh];
        else
            R.peakamp(Nspect,1) = (SpectAllF(1,2)-SpectAllF(1,1))*Spect(1,T(Nspect,1).peakfreq(1));
            T.peakamp(Nspect,:) = T(Nspect,1).peakfreq;
            R.peakratio(Nspect,1) = (SpectAllF(1,2)-SpectAllF(1,1))*Spect(1,T(Nspect,1).peakfreq(1)) / ...
                trapz(SpectAllF(1,freqlow:freqhigh),Spect2,2);
            T.peakratio(Nspect,:) = [freqlow freqhigh];
            
        end
        if isfield(settings,'freqs') & ~isempty(settings.freqs)
            tmp = find(SpectAllF >= min(settings.freqs(:,1)),1);
        else
            tmp = 1;
        end
        [minfreq,minfreqI] = min(Spect(tmp:T.peakfreq(Nspect,1)));
        minfreqI = minfreqI + tmp - 1;
        R.peak2min(Nspect,1) = Spect(T.peakfreq(Nspect,1)) / Spect(minfreqI(1));
        T.peak2min(Nspect,:) = [minfreqI minfreq];
        clearvars minfreq minfreqI tmp
    end
    
    % calc median freq
    freqlow = find(SpectAllF >= settings.lowfreqcog, 1 );
    freqhigh = find(SpectAllF <= settings.highfreqcog, 1,'last');
    if ~isempty(freqlow) & ~isempty(freqhigh)
        Spect2 = Spect(1,freqlow:freqhigh);
        [rc,cc] = ndgrid(1:size(Spect2,1),1:size(Spect2,2));
        Spectt = sum(Spect2(:));
        c1 = sum(Spect2(:) .* rc(:)) / Spectt;
        c2 = sum(Spect2(:) .* cc(:)) / Spectt;
        R.cogfreqAll(Nspect,1) = SpectAllF(freqlow -1) + (c2 * (SpectAllF(2) - SpectAllF(1)));
        tmp = ((R.cogfreq(Nspect,1) - SpectAllF(1))/(SpectAllF(2)-SpectAllF(1)))+1;
        T.cogfreq(Nspect,:) = [tmp tmp];
        R.cogareapower(Nspect,1) = trapz(SpectAllF(1,freqlow:freqhigh),Spect2,2);
        T.cogareapower(Nspect,:) = [freqlow freqhigh];
    end
    clearvars Spect2 Spectt maxindex maxvalue c1 c2 rc cc
    
    % Calculate Power last frequency band
    if isfield(settings,'freqs') & ~isempty(settings.freqs) & doquality == true
        freqlow = find(SpectAllF >= settings.freqs(end,1), 1 );
        freqhigh = find(SpectAllF <= settings.freqs(end,2), 1,'last');
        R.bplast(Nspect,1) = trapz(SpectAllF(1,freqlow:freqhigh),Spect(1,freqlow:freqhigh),2);
        T.bplast(Nspect,:) = [freqlow freqhigh];
    end
end

% Calculate Median Frequency
[cogfreq,Tcogfreq] = lab_collect_spectraldata_mf(SpectAllM,SpectAllF,settings);
R.cogfreq = cogfreq;
T.cogfreq = Tcogfreq;
clearvars cogfreq Tcogfreq

if settings.domedian == true
    T.cogareapower = median(T.cogareapower,1);
    T.areapower = median(T.areapower,1);
    T.bplast = median(T.bplast,1);
    R.cogareapower = median(R.cogareapower,1);
    R.areapower = median(R.areapower,1);
    R.bplast = median(R.bplast,1);
    if sum(isnan(R.peakfreqAll)) <= length(R.peakfreqAll)/2
        T.peakfreq = median(T.peakfreq(~isnan(R.peakfreqAll),:),1);
        T.peakamp = median(T.peakamp(~isnan(R.peakamp),:),1);
        T.peakratio = median(T.peakratio(~isnan(R.peakratio),:),1);
        T.peak2min = median(T.peak2min(~isnan(R.peak2min),:),1);
        R.peakfreq = median(R.peakfreqAll(~isnan(R.peakfreqAll)),1);
        R.peakamp = median(R.peakamp(~isnan(R.peakamp)),1);
        R.peakratio = median(R.peakratio(~isnan(R.peakratio)),1);
        R.peak2min = median(R.peak2min(~isnan(R.peak2min)),1);
    else
        R.peakfreq = NaN;
        R.peakamp = 0;
        R.peakratio = 0;
        R.peak2min = 0;
        T.peakfreq = zeros(1,2);
        T.peakamp = zeros(1,2);
        T.peakratio = zeros(1,2);
        T.peak2min = zeros(1,2);
    end
    T.Spect = median(SpectAllM,1);
else
    T.cogareapower = mean(T.cogareapower,1);
    T.areapower = mean(T.areapower,1);
    T.bplast = mean(T.bplast,1);
    R.cogareapower = mean(R.cogareapower,1);
    R.areapower = mean(R.areapower,1);
    R.bplast = mean(R.bplast,1);
    if sum(isnan(R.peakfreqAll)) <= length(R.peakfreqAll)/2
        T.peakfreq = mean(T.peakfreq(~isnan(R.peakfreqAll),:),1);
        T.peakamp = mean(T.peakamp(~isnan(R.peakamp),:),1);
        T.peakratio = mean(T.peakratio(~isnan(R.peakratio),:),1);
        T.peak2min = mean(T.peak2min(~isnan(R.peak2min),:),1);
        R.peakfreq = mean(R.peakfreqAll(~isnan(R.peakfreqAll)),1);
        R.peakamp = mean(R.peakamp(~isnan(R.peakamp)),1);
        R.peakratio = mean(R.peakratio(~isnan(R.peakratio)),1);
        R.peak2min = mean(R.peak2min(~isnan(R.peak2min)),1);
    else
        R.peakfreq = NaN;
        R.peakamp = 0;
        R.peakratio = 0;
        R.peak2min = 0;
        T.peakfreq = zeros(1,2);
        T.peakamp = zeros(1,2);
        T.peakratio = zeros(1,2);
        T.peak2min = zeros(1,2);
    end
    T.Spect = mean(SpectAllM .^0.5,1).^2;
end
if isempty(R.peakfreq)
    R.peakfreq = NaN;
    T.peakfreq = [0 0];
    R.peakamp = NaN;
    T.peakamp = [0 0];
    R.peakratio = NaN;
    T.peakratio = [0 0];
    R.peak2min = NaN;
    T.peak2min = [0 0];
end

if ~isempty(header)
    if isfield(header,'numchannels')
        [T,settings] = lab_test_badchannels(header,settings,T);
    else
        if isfield(header,'badchans')
            T.badchans = header.badchans;
        else
            T.badchans = 0;
        end
        if isfield(header,'badact')
            T.badact = header.badact;
        else
            T.badact = 0;
        end
    end
else
    T.badchans = 0;
    T.badact = 0;
end
R.badchans = T.badchans;
R.badact = T.badact;

if isfield(header,'quality')
    T.epochquality = mean(header.quality);
else
    T.epochquality = 0;
end
R.epochquality = T.epochquality;

if isfield(settings,'QUALITY') & ~isempty(settings.QUALITY)
    R = lab_evaluate_quality(R,settings.QUALITY);
end
if ~isfield(R,'bchans')
    R.bchans = sum(R.badchans>0);
end
if ~isfield(R,'good')
    R.good = false;
end

end

function [PF,maxindex] = calcPF(Spect,SpectAllF,freqlow,freqhigh,MinPeak2Min)

numchans = size(Spect,1);
PF = NaN(numchans,1);
maxindex = zeros(numchans,1);

if isempty(freqlow) | isempty(freqhigh)
    return
end

for Nspect = 1:numchans
    P2M = 0;
    Lag = 0;
    MaxLag = freqhigh - freqlow;
    while P2M < MinPeak2Min & Lag < MaxLag
        [maxvalue,maxindex] = max(Spect(Nspect,freqlow+Lag:freqhigh),[],2);
        maxindex = maxindex + Lag;
        if length(maxindex) > 1
            tmp = find(maxindex>find(SpectAllF >= 8,1));
            if ~isempty(tmp)
                maxindex = maxindex(tmp(1));
            else
                maxindex = maxindex(end);
            end
            clearvars tmp
        end
        minvalue = min(Spect(Nspect,1:maxindex),[],2);
        P2Mtmp = maxvalue/minvalue;
        Lag = Lag + 1;
        if P2Mtmp >= P2M
            P2M = P2Mtmp;
            
        end
    end
    clearvars P2M Lag MaxLag
    
    if ~isempty(maxvalue) & maxvalue > 0 & maxindex > 1 & maxindex < size(Spect,2)
        maxindex = (freqlow - 1 + round(maxindex));
        PF(Nspect,1) = SpectAllF(maxindex);
    end
end

end
