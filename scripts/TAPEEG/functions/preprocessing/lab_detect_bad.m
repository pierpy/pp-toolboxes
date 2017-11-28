% Detects bad channels: combination of txt-input, FASTER, Fieldtrip and TAPEEG routines
%
% [badelectrodes,bad,settings] = lab_detect_bad(data,header,settings,cfg)
%
% data     = eeg/meg data (chans x timeframes)
% header   = output of lab_read_data
% settings = structure with config (optional)
%
%
% Written by F. Hatz 2014

function [badelectrodes,bad,settings] = lab_detect_bad(data,header,settings,cfg,noverbose,silent)

if ~exist('data','var')
    [data,header,cfg] = lab_read_data;
    if isempty(data)
        return
    end
end
if ~exist('header','var') | ~isfield(header,'samplingrate')
    header = lab_create_header(data);
end
if ~exist('cfg','var')
    cfg = [];
end
if ~exist('settings','var')
    settings = [];
    [settings,skipprocessing] = lab_set_detect_bad(settings,header,cfg,0,0,2,0,1,1,1,1,1);
    if skipprocessing == 1
        badelectrodes = [];
        bad = [];
        return
    end
end
if ~isfield(settings,'length') | isempty(settings.length)
    settings.length = 4;
end
wlength = ceil(settings.length * header.samplingrate);
if wlength > size(data,2)
    wlength = size(data,2);
end

if isfield(header,'highpass') & max(settings.spectslow) < header.highpass
    settings.spectslow = [0 0];
end
if isfield(header,'lowpass') & min(settings.spectshigh) > header.lowpass
    settings.spectshigh = [0 0];
end

if isfield(settings,'filemethod')
    [bad,settings,skipautomated] = lab_get_filebad(header,settings,cfg);
else
    skipautomated = 0;
    bad = [];
end
if isfield(bad,'file')
    badelectrodes = bad.file;
else
    badelectrodes = [];
end

if isfield(settings,'fixedbad') & ~isempty(settings.fixedbad)
    badelectrodes = union(badelectrodes,settings.fixedbad);
    bad.fixed = settings.fixedbad(:)';
end

if skipautomated == 1
    return
end

bad.error = 1;
if ~isfield(settings,'percentbad')
    percentbad = 1;
elseif settings.percentbad > 100
    bad.error = 0;
    disp('   Automated bad channels search disabled (percent bad > 100)')
    return
else
    percentbad = settings.percentbad / 100;
end

% Prepare data
if ~isfield(settings,'blinkchans') & isfield(settings,'eog') & ~isempty(settings.eog)
    header.bad.All = badelectrodes;
    settings = lab_calculate_eog(data,header,settings);
end
if isfield(header,'numdatachannels') & size(data,1) > header.numdatachannels
    data = data(1:header.numdatachannels,:);
end
data = data(:,1:wlength * floor(size(data,2)/wlength));
data = reshape(data,size(data,1),wlength,floor(size(data,2)/wlength));
datacounter = floor(size(data,3)/20);
if datacounter < 1
    datacounter = 1;
end

if isfield(settings,'freqlim60') & ~isempty(settings.freqlim60) & settings.freqlim60 > 0
    dospect = 1;
elseif isfield(settings,'freqlim50') & ~isempty(settings.freqlim50) & settings.freqlim50 > 0
    dospect = 1;
elseif isfield(settings,'freqlimlow') & ~isempty(settings.freqlimlow) & settings.freqlimlow(1,1) > 0
    dospect = 1;
elseif isfield(settings,'freqlimhigh') & ~isempty(settings.freqlimhigh) & settings.freqlimhigh(1,1) > 0
    dospect = 1;
else
    dospect = 0;
end

if dospect == 1
    % Calculate spectras
    if ~exist('silent','var')
        fprintf('     calculate spectras for bad channels analysis')
    end
    hanning = hann(ceil(size(data,2) * 0.2))';
    window = ones(1,ceil(size(data,2)));
    window(1:ceil(size(hanning,2)/2)) = hanning(1:ceil(size(hanning,2)/2));
    window(end-ceil(size(hanning,2)/2):end) = hanning(end-ceil(size(hanning,2)/2):end);
    window = repmat(window,size(data,1),1)';
    for indextrials = 1:size(data,3)
        [specdata(:,:,indextrials)] = abs(fft(detrend(data(:,:,indextrials)').*window,size(window,1),1)).^2'; %#ok<AGROW>
        if mod(indextrials,datacounter) == 0 & ~exist('silent','var')
            fprintf('.')
        end
    end
    specdata = specdata(:,1:ceil(size(window,1)/2)+1,:);
    freqs = header.samplingrate/2*linspace(0,1,ceil(size(window,1)/2)+1);
    for i = 1:size(specdata,3)
        specdata(:,:,i) = specdata(:,:,i) ./ repmat(sum(specdata(:,1:find(freqs > 40, 1 ),i),2),1,size(specdata,2));
    end
    clearvars window
    if ~exist('silent','var')
        disp(':')
    end
end

bad.epochs = zeros(size(data,1),size(data,3));
if ~isempty(badelectrodes)
    bad.epochs(badelectrodes,:) = 1;
    data(badelectrodes,:,:) = NaN;
end

% Find broken channels
if isfield(settings,'zvaluebroken') & settings.zvaluebroken > 0
    if ~exist('silent','var')
        disp('     detect broken channels by absolute power')
    end
    badbroken = zeros(size(data,1),size(data,3));
    for i = 1:size(data,3)
        SumPow = sum(abs(data(:,:,i)),2);
        flagtmp = 1;
        while ~isempty(flagtmp)
            tmp = (SumPow - median(SumPow)) / std(SumPow);
            flagtmp = find(tmp > settings.zvaluebroken);
            if ~isempty(flagtmp)
                badbroken(flagtmp,i) = 1;
                bad.epochs(flagtmp,i) = 1;
                SumPow(flagtmp) = median(SumPow);
            end
        end
    end
    clearvars tmp flagtmp
    rejE = sum(badbroken,2) ./ size(badbroken,2);
    rejE = find(rejE > percentbad)';
    SumPow = sum(abs(reshape(data,[size(data,1) size(data,2)*size(data,3)])),2);
    flagtmp = 1;
    badbroken = zeros(size(data,1),1);
    while ~isempty(flagtmp)
        tmp = (SumPow - median(SumPow)) / std(SumPow);
        flagtmp = find(tmp > settings.zvaluebroken);
        if ~isempty(flagtmp)
            badbroken(flagtmp,1) = 1;
            SumPow(flagtmp) = median(SumPow);
        end
    end
    rejE = union(rejE,find(badbroken == 1));
    badelectrodes = union(badelectrodes,rejE);
    bad.broken = rejE;
    clearvars rejE badbroken
    if ~isempty(badelectrodes)
        data(badelectrodes,:,:) = NaN;
    end
end

% Look for spectral artifact 50Hz
if isfield(settings,'freqlim50') & ~isempty(settings.freqlim50) & settings.freqlim50 > 0
    bad.epochsSpect50 = zeros(size(data,1),size(data,3));
    bad.raw.epochsSpect50 = zeros(size(data,1),size(data,3));
    badepochs = zeros(size(specdata,1),size(specdata,3));
    for i = 1:size(specdata,1)
        for j = 1:size(specdata,3)
            bad.raw.epochsSpect50(i,j) = sum(specdata(i,find(freqs > 49,1):find(freqs > 51,1),j));
            if bad.raw.epochsSpect50(i,j) > settings.freqlim50/100
                bad.epochs(i,j) = 1;
                bad.epochsSpect50(i,j) = 1;
                badepochs(i,j) = 1;
            end
        end
    end
    rejE = sum(badepochs,2) ./ size(badepochs,2);
    rejE = find(rejE > percentbad)';
    badelectrodes = union(badelectrodes,rejE);
    bad.spectral50 = rejE;
    clearvars rejE badepochs
else
    bad.spectral50 = [];
end

% Look for spectral artifact 60Hz
if isfield(settings,'freqlim60') & ~isempty(settings.freqlim60) & settings.freqlim60 > 0
    bad.epochsSpect60 = zeros(size(data,1),size(data,3));
    bad.raw.epochsSpect60 = zeros(size(data,1),size(data,3));
    badepochs = zeros(size(specdata,1),size(specdata,3));
    for i = 1:size(specdata,1)
        for j = 1:size(specdata,3)
            bad.raw.epochsSpect60(i,j) = sum(specdata(i,find(freqs > 59,1):find(freqs > 61,1),j));
            if bad.raw.epochsSpect60(i,j) > settings.freqlim60/100
                bad.epochs(i,j) = 1;
                bad.epochsSpect60(i,j) = 1;
                badepochs(i,j) = 1;
            end
        end
    end
    rejE = sum(badepochs,2) ./ size(badepochs,2);
    rejE = find(rejE >= percentbad)';
    badelectrodes = union(badelectrodes,rejE);
    bad.spectral60 = rejE;
    clearvars rejE badepochs
else
    bad.spectral60 = [];
end

% Look for spectral artifact highfreq
if isfield(settings,'freqlimhigh') & ~isempty(settings.freqlimhigh) & settings.freqlimhigh > 0 &settings.spectshigh(1,1) > 0
    bad.epochsSpectHigh = zeros(size(data,1),size(data,3));
    bad.raw.epochsSpectHigh = zeros(size(data,1),size(data,3));
    badepochs = zeros(size(specdata,1),size(specdata,3));
    for i = 1:size(specdata,1)
        for j = 1:size(specdata,3)
            bad.raw.epochsSpectHigh(i,j) = sum(specdata(i,find(freqs > settings.spectshigh(1,1),1):find(freqs > settings.spectshigh(1,2),1),j));
            if bad.raw.epochsSpectHigh(i,j) > settings.freqlimhigh/100
                bad.epochs(i,j) = 1;
                bad.epochsSpectHigh(i,j) = 1;
                badepochs(i,j) = 1;
            end
        end
    end
    rejE = sum(badepochs,2) ./ size(badepochs,2);
    rejE = find(rejE >= percentbad)';
    badelectrodes = union(badelectrodes,rejE);
    bad.spectralhigh = rejE;
    clearvars rejE badepochs
else
    bad.spectralhigh = [];
end

% Look for spectral artifact lowfreq
if isfield(settings,'freqlimlow') & ~isempty(settings.freqlimlow) & settings.freqlimlow > 0 & settings.spectslow(1,2) > 0
    bad.epochsSpectLow = zeros(size(data,1),size(data,3));
    bad.raw.epochsSpectLow = zeros(size(data,1),size(data,3));
    badepochs = zeros(size(specdata,1),size(specdata,3));
    for i = 1:size(specdata,1)
        bad.raw.epochsSpectLow(i,j) = sum(specdata(i,find(freqs > settings.spectslow(1,1),1):find(freqs > settings.spectslow(1,2),1),j));
        for j = 1:size(specdata,3)
            if bad.raw.epochsSpectLow(i,j) > settings.freqlimlow/100
                bad.epochs(i,j) = 1;
                bad.epochsSpectLow(i,j) = 1;
                badepochs(i,j) = 1;
            end
        end
    end
    rejE = sum(badepochs,2) ./ size(badepochs,2);
    rejE = find(rejE >= percentbad)';
    badelectrodes = union(badelectrodes,rejE);
    bad.spectrallow = rejE;
    clearvars rejE badepochs
else
    bad.spectrallow = [];
end

% Calculate variance of channels
if isfield(settings,'zvaluevars') & ~isempty(settings.zvaluevars) & settings.zvaluevars > 0
    bad.epochsVar = zeros(size(data,1),size(data,3));
    if ~exist('silent','var')
        disp('     calculate variance for bad channels analysis')
    end
    % Read electrode localisations
    eeg_chans = 1:size(data,1);
    if isfield(header,'ref_chan') & ~isempty(header.ref_chan) & isnumeric(header.ref_chan) & ...
            length(header.ref_chan)==1 & header.ref_chan ~= 0 & isfield(header,'locs') & ~isempty(header.locs)
        pol_dist = lab_distancematrix(data,header.locs);
        [s_pol_dist,dist_inds] = sort(pol_dist(header.ref_chan,eeg_chans));
        [~, idist_inds] = sort(dist_inds);
    end
    for i = 1:size(data,3)
        vars = var(data(:,:,i),[],2)';
        vars(~isfinite(vars)) = mean(vars(isfinite(vars)));
        if exist('s_pol_dist','var')
            p = polyfit(s_pol_dist,vars(dist_inds),2);
            fitcurve = polyval(p,s_pol_dist);
            corrected = vars - fitcurve(idist_inds);
            rejE(:,i) = corrected;
        else
            rejE(:,i) = vars;
        end
    end
    for u = 1:size(rejE,2)
        rejE(isnan(rejE(:,u)),u) = nanmean(rejE(:,u));
        rejE(:,u) = rejE(:,u) - median(rejE(:,u));
    end
    zs = rejE - repmat(mean(rejE,1),size(rejE,1),1);
    zs = zs./repmat(std(zs,[],1),size(rejE,1),1);
    zs(isnan(zs))=0;
    bad.raw.epochsVar = abs(zs);
    rejE = abs(zs) > settings.zvaluevars;
    bad.epochs(rejE == 1) = 1;
    bad.epochsVar(rejE == 1) = 1;
    rejE = sum(rejE,2) ./ size(rejE,2);
    bad.variance = find(rejE >= percentbad)';
    badelectrodes = union(badelectrodes,bad.variance);
    clearvars rejE zs
else
    bad.variance = [];
end

% Calculate hurst exponent
if isfield(settings,'zvaluehurst') & ~isempty(settings.zvaluehurst) & settings.zvaluehurst > 0
    bad.epochsHurst = zeros(size(data,1),size(data,3));
    if ~exist('silent','var')
        fprintf('     calculate hurst exponent for bad channels analysis')
    end
    for i = 1:size(data,3)
        for j = 1:size(data,1)
            rejE(j,i) = hurst_exponent(data(j,:,i));
        end
        if mod(i,datacounter) == 0 & ~exist('silent','var')
            fprintf('.')
        end
    end
    for u = 1:size(rejE,2)
        rejE(isnan(rejE(:,u)),u)=nanmean(rejE(:,u));
        rejE(:,u) = rejE(:,u) - median(rejE(:,u));
    end
    zs = rejE - repmat(mean(rejE,1),size(rejE,1),1);
    zs = zs./repmat(std(zs,[],1),size(rejE,1),1);
    zs(isnan(zs))=0;
    bad.raw.epochsHurst = abs(zs);
    rejE = abs(zs) > settings.zvaluehurst;
    bad.epochs(rejE == 1) = 1;
    bad.epochsHurst(rejE == 1) = 1;
    rejE = sum(rejE,2) ./ size(rejE,2);
    bad.hurst = find(rejE >= percentbad)';
    badelectrodes = union(badelectrodes,bad.hurst);
    clearvars rejE zs
    if ~exist('silent','var')
        disp(':')
    end
else
    bad.hurst = [];
end

% Check for markers to exclude
if isfield(header,'events') & isfield(header.events,'POS') & ~isempty(header.events.POS) & isfield(settings,'markerexclude') & ~isempty(settings.markerexclude)
    if ~exist('silent','var')
        disp('     define bad-channel-epochs by markers')
    end
    bad.epochsMarker = zeros(size(data,1),size(data,3));
    markerdata = zeros(1,size(data,2)*size(data,3));
    if isempty(settings.markerexclude)
        markers = [];
    elseif strcmp(settings.markerexclude{1,1},'all')
        markers = unique(header.events.TYP);
    else
        markers = settings.markerexclude;
    end
    for i = 1:length(markers)
        tmp = find(strcmp(header.events.TYP,markers(i)));
        if ~isempty(tmp)
            for j = 1:size(tmp,2)
                markerdata(1,header.events.POS(tmp(1,j)):header.events.POS(tmp(1,j))+header.events.DUR(tmp(1,j))-1) = 1;
            end
        end
    end
    clearvars tmp
    markerdata = markerdata(1,1:size(data,2)*size(data,3));
    markerdata = reshape(markerdata,size(data,2),size(data,3));
    markerdata = sum(markerdata,1) > 0;
    bad.epochs(:,markerdata) = 1;
    bad.epochsMarker(:,markerdata) = 1;
end

% find corr to eye channels
if isfield(settings,'zvalueeye') & ~isempty(settings.zvalueeye) & settings.zvalueeye > 0
    if isfield(settings,'blinkchans') & ~isempty(settings.blinkchans)
        bad.epochsEye = zeros(size(data,1),size(data,3));
        if ~exist('silent','var')
            disp('     detect correlations to eye channels')
        end
        [B, A] = butter(4,[2/header.samplingrate 30/header.samplingrate]);
        for v=1:size(settings.blinkchans,1)
            blink_chans(v,:)=filtfilt(B,A,settings.blinkchans(v,:)); %#ok<AGROW>
        end
        blink_chans = blink_chans(:,1:wlength * floor(size(blink_chans,2)/wlength));
        blink_chans = reshape(blink_chans,size(blink_chans,1),wlength,floor(size(blink_chans,2)/wlength));
        for i = 1:size(data,3)
            for j = 1:size(data,1)
                tmp = filtfilt(B,A,data(j,:,i));
                for v = 1:size(blink_chans,1)
                    if max(blink_chans(v,:,i))~=0 & min(blink_chans(v,:,i))~=0;
                        f = corrcoef(tmp,blink_chans(v,:,i));
                        x(v) = abs(f(1,2)); %#ok<AGROW>
                    else
                        x(v) = v; %#ok<AGROW>
                    end
                end
                clearvars tmp
                rejE(j,i) = max(x);
            end
        end
        zs = rejE - repmat(mean(rejE,1),size(rejE,1),1);
        zs=zs./repmat(std(zs,[],1),size(rejE,1),1);
        zs(isnan(zs))=0;
        bad.raw.epochsEye = abs(zs);
        rejE = abs(zs) > settings.zvalueeye;
        bad.epochs(rejE == 1) = 1;
        bad.epochsEye(rejE == 1) = 1;
        rejE = sum(rejE,2) ./ size(rejE,2);
        bad.eye = find(rejE >= percentbad)';
        badelectrodes = union(badelectrodes,bad.eye);
        clearvars rejE zs x
    else
        if ~exist('silent','var')
            disp('     skip eye, no eog defined')
        end
        bad.eye = [];
    end
else
    bad.eye = [];
end

% Calculate median
if isfield(settings,'zvaluemedian') & ~isempty(settings.zvaluemedian) & settings.zvaluemedian > 0
    bad.epochsMedian = zeros(size(data,1),size(data,3));
    if ~exist('silent','var')
        disp('     calculate median gradient for bad channels analysis')
    end
    for i = 1:size(data,3)
        for j = 1:size(data,1)
            rejE(j,i) = median(diff(data(j,:,i)));
        end
    end
    for u = 1:size(rejE,2)
        rejE(isnan(rejE(:,u)),u)=nanmean(rejE(:,u));
        rejE(:,u) = rejE(:,u) - median(rejE(:,u));
    end
    zs = rejE - repmat(mean(rejE,1),size(rejE,1),1);
    zs=zs./repmat(std(zs,[],1),size(rejE,1),1);
    zs(isnan(zs))=0;
    bad.raw.epochsMedian = abs(zs);
    rejE = abs(zs) > settings.zvaluemedian;
    bad.epochs(rejE == 1) = 1;
    bad.epochsMedian(rejE == 1) = 1;
    rejE = sum(rejE,2) ./ size(rejE,2);
    bad.median = find(rejE >= percentbad)';
    badelectrodes = union(badelectrodes,bad.median);
    clearvars rejE zs
else
    bad.median = [];
end

% Calculate amplitude
if isfield(settings,'zvalueamplitude') & ~isempty(settings.zvalueamplitude) & settings.zvalueamplitude > 0
    bad.epochsAmplitude = zeros(size(data,1),size(data,3));
    if ~exist('silent','var')
        disp('     calculate amplitude gradient for bad channels analysis')
    end
    for i = 1:size(data,3)
        for j = 1:size(data,1)
            rejE(j,i) = max(data(j,:,i)) - min(data(j,:,i));
        end
    end
    for u = 1:size(rejE,2)
        rejE(isnan(rejE(:,u)),u)=nanmean(rejE(:,u));
        rejE(:,u) = rejE(:,u) - median(rejE(:,u));
    end
    zs = rejE - repmat(mean(rejE,1),size(rejE,1),1);
    zs = zs./repmat(std(zs,[],1),size(rejE,1),1);
    zs(isnan(zs))=0;
    bad.raw.epochsAmplitude = abs(zs);
    rejE = abs(zs) > settings.zvalueamplitude;
    bad.epochs(rejE == 1) = 1;
    bad.epochsAmplitude(rejE == 1) = 1;
    rejE = sum(rejE,2) ./ size(rejE,2);
    bad.amplitude = find(rejE >= percentbad)';
    badelectrodes = union(badelectrodes,bad.amplitude);
    clearvars rejE zs
else
    bad.amplitude = [];
end

% Calculate kurtosis
if isfield(settings,'zvaluekurtosis') & ~isempty(settings.zvaluekurtosis) & settings.zvaluekurtosis > 0
    if isfield(header,'W');
        bad.epochsKurtosis = zeros(size(data,1),size(data,3));
        if ~exist('silent','var')
            disp('     calculate kurtosis for bad channels analysis')
        end
        for i = 1:size(header.W,2)
            rejE(i,1) = kurt(header.W(:,i));
        end
        zs = rejE - repmat(mean(rejE,1),size(rejE,1),1);    
        zs=zs./repmat(std(zs,[],1),size(rejE,1),1);
        zs(isnan(zs))=0;
        bad.raw.epochsKurtosis = abs(zs);
        rejE = abs(zs) > settings.zvaluekurtosis;
        bad.epochs(rejE == 1,:) = 1;
        bad.epochsKurtosis(rejE == 1,:) = 1;
        bad.kurtosis = find(rejE == 1)';
        badelectrodes = union(badelectrodes,bad.kurtosis);
        clearvars rejE zs
    elseif ~exist('silent','var')
        disp('     skip kurtosis, no information on weights')
    end
else
    bad.kurtosis = [];
end

% Calculate topo correlation
if isfield(settings,'zvaluecorr') & ~isempty(settings.zvaluecorr) & settings.zvaluecorr > 0
    if isfield(header,'locs') & ~isempty(header.locs) & length(header.locs.x) == size(data,1);
        bad.epochsCorr = zeros(size(data,1),size(data,3));
        if ~exist('silent','var')
            disp('     calculate topo correlation for bad channels analysis')
        end
        if ~isfield(settings,'LAPL')
            settings.LAPL = [];
        end
        [refmatrix,~,settings.LAPL] = lab_calc_lapmatrix(header,settings.LAPL);
        refmatrix(1:size(refmatrix,1)+1:end) = 0;
        for i = 1:size(data,3)
            reftmp = refmatrix * data(:,:,i);
            for j = 1:size(data,1);
                rejE(j,i) = abs(corr(data(j,:,i)',reftmp(j,:)'));
            end
        end
        for u = 1:size(rejE,2)
            rejE(isnan(rejE(:,u)),u) = nanmean(rejE(:,u));
            rejE(:,u) = rejE(:,u) - median(rejE(:,u));
        end
        zs = rejE - repmat(mean(rejE,1),size(rejE,1),1);
        zs = zs./repmat(std(zs,[],1),size(rejE,1),1);
        zs(isnan(zs))=0;
        bad.raw.epochsCorr = abs(zs);
        rejE = abs(zs) > settings.zvaluecorr;
        bad.epochs(rejE == 1) = 1;
        bad.epochsCorr(rejE == 1) = 1;
        rejE = sum(rejE,2) ./ size(rejE,2);
        bad.corr = find(rejE >= percentbad)';
        badelectrodes = union(badelectrodes,bad.corr);
        clearvars rejE zs
    else
        if ~exist('silent','var')
            disp('     skip correlation, no electrodes-locs')
        end
    end
else
    bad.corr = [];
end

% detect minimal peak2min
if isfield(settings,'PEAK2MIN') & isfield(settings.PEAK2MIN,'BAchannels') & ~isempty(settings.PEAK2MIN.BAchannels)
    if isfield(header,'includechans') & ~isempty(header.includechans)
        tmp = zeros(1,max(max(settings.PEAK2MIN.BAchannels),max(header.includechans)));
        tmp(settings.PEAK2MIN.BAchannels) = 1;
        tmp = tmp(header.includechans);
        BAchannels = find(tmp==1);
        clearvars tmp
    elseif isfield(header,'montage') & header.montage == true
        BAchannels = 1:size(data,1);
    else
        BAchannels = settings.PEAK2MIN.BAchannels;
    end
    if max(BAchannels) > size(data,1)
        BAchannels = 1:size(data,1);
    end
    fprintf('     calculate peak2min');
    peak2min = ones(1,size(data,3));
    for indextrials = 1:size(data,3)
        [spectra,~,freqs] = pmtm(detrend(nanmean(data(BAchannels,:,indextrials),1)),[],[],header.samplingrate);
        spectra = spectra';
        freqs = freqs';
        if mod(indextrials,datacounter) == 0 & ~exist('silent','var')
            fprintf('.')
        end
        if ~exist('freqlow','var')
            freqlow = find(freqs >= settings.PEAK2MIN.lowfreqpeak, 1 );
            freqhigh = find(freqs <= settings.PEAK2MIN.highfreqpeak, 1,'last');
        end
        if ~isempty(freqlow) & ~isempty(freqhigh)
            spectra2 = spectra(1,freqlow:freqhigh);
            P2M = 0;
            Lag = 0;
            MaxLag = size(spectra2,2);
            while P2M < settings.PEAK2MIN.MinPeak2Min & Lag < MaxLag
                [maxvalue,maxindex] = max(spectra2(1,1+Lag:end),[],2);
                maxindex = maxindex + Lag;
                if length(maxindex) > 1
                    tmp = find(maxindex>find(freqs >= 8,1));
                    if ~isempty(tmp)
                        maxindex = maxindex(tmp(1));
                    else
                        maxindex = maxindex(end);
                    end
                    clearvars tmp
                end
                minvalue = min(spectra2(1:maxindex),[],2);
                P2M = maxvalue/minvalue;
                Lag = Lag + 1;
                if P2M >= settings.PEAK2MIN.MinPeak2Min | peak2min(1,indextrials) == 1
                    peak2min(1,indextrials) = P2M;
                end
            end
            clearvars Lag MaxLag P2M maxvalue maxindex minvalue spectra2
        end
        % Spect(indextrials,:) = spectra;
        clearvars spectra
    end
    if ~exist('silent','var')
        disp(':');
    end
    clearvars freqlow freqhigh
    bad.raw.peak2min = peak2min;
    if strcmp(settings.PEAK2MIN.mode,'threshold')
        if isempty(settings.PEAK2MIN.threshold)
            settings.PEAK2MIN.threshold = 1.5;
        end
        if isempty(settings.PEAK2MIN.factor)
            settings.PEAK2MIN.threshold = 2;
        end
        bad.peak2min = ones(1,size(data,3));
        bad.peak2min(peak2min < settings.PEAK2MIN.threshold) = 1/settings.PEAK2MIN.factor;
    else
        if isempty(settings.PEAK2MIN.factor)
            settings.PEAK2MIN.threshold = 3;
        end
        peak2min = peak2min / ((max(peak2min) - min(peak2min)) / (1 - 1/settings.PEAK2MIN.factor));
        bad.peak2min = peak2min - max(peak2min) + 1;
    end
else
    bad.peak2min = [];
end

% evaluate microstates correlation
if isfield(settings,'MicroCorr') & settings.MicroCorr == true & isfield(header,'CORR') & ~isempty(header.CORR)
    disp('     get microstates correlations')
    CORR = reshape(header.CORR(1,1:size(data,2)*size(data,3)),size(data,2),size(data,3));
    bad.microcorr = mean(abs(CORR),1);
    if isfield(settings,'MicroCorrNeg') & settings.MicroCorrNeg == true
        disp('       convert to negativ microstates correlations')
        bad.microcorr = 1 - bad.microcorr;
    end
else
    bad.microcorr = [];
end

% do analysis without epochs
bad.window = size(data,2);
data = reshape(data,size(data,1),size(data,2)*size(data,3));

% find channels with ecg artifact
if isfield(settings,'ecgdetect') & ~isempty(settings.ecgdetect) & settings.ecgdetect < 3
    if ~exist('silent','var')
        disp('     detect activations with ecg artifact')
    end
    if settings.ecgdetect == 1
        if isfield(settings,'ecg_chan') & length(settings.ecg_chan) >= size(data,2)
            settings.ecg_chan = settings.ecg_chan(1,1:size(data,2));
        elseif isfield(header,'ecg_ch') & header.ecg_ch > 0 & header.ecg_ch <= size(data,1)
            settings.ecg_chan = data(header.ecg_ch,:);
        elseif isfield(settings,'ecg_ch') & settings.ecg_ch > 0 & settings.ecg_ch <= size(data,1)
            settings.ecg_chan = data(settings.ecg_ch,:);
        else
            settings.ecg_chan = [];
            if ~exist('silent','var')
                disp('      detect ecg disabled, wrong ecg-channel information')
            end
        end
        bad.QRS=[];
        if ~isempty(settings.ecg_chan)
            ECGchloc =lab_nqrsdetect(settings.ecg_chan,header.samplingrate)';
            if length(ECGchloc) > size(data,2)/(header.samplingrate*1.5)
                ECGch = zeros(1,size(data,2));
                for cp = -round(header.samplingrate/100):round(header.samplingrate/10)
                    ECGchloc2 = ECGchloc + cp;
                    ECGch(1,ECGchloc2) = 1;
                end
                clearvars ECGchloc ECGchloc2;
                for cp = 1:size(data,1);
                    tmp =lab_nqrsdetect(data(cp,:),header.samplingrate)';
                    if sum(ECGch(1,tmp)) > length(tmp)*percentbad
                        bad.QRS = [bad.QRS cp];
                    end
                end
                clearvars tmp ECGch;
                badelectrodes = union(badelectrodes,bad.QRS);
            end
        end
    elseif settings.ecgdetect == 2
        if ~isfield(settings,'ecg_ch') | length(settings.ecg_ch) ~= 2
            settings.ecg_ch = [50 120];
            if ~exist('silent','var')
                disp('      wrong ecg-rate information, set to 50-120')
            end
        end
        for cp = 1:size(data,1)
            bad.QRS(cp) = 0;
            tmp =lab_nqrsdetect(data(cp,:),header.samplingrate)';
            if length(tmp) > 20
                tmp = tmp(2:end) - tmp(1:end-1);
                tmp = sort(tmp);
                tmp = tmp(ceil(0.1*length(tmp)):floor(0.9*length(tmp)));
                if (mean(tmp) > (header.samplingrate*(60/settings.ecg_ch(1,2)))) & ...
                        (mean(tmp) < (header.samplingrate*(60/settings.ecg_ch(1,1)))) & ...
                        (std(tmp) < (header.samplingrate/8)) & (size(tmp,2) > (size(data,2)/(mean(tmp)*4)))
                    bad.QRS(cp) = 1;
                end
            end
            clearvars tmp
        end
        bad.QRS = find(bad.QRS);
        badelectrodes = union(badelectrodes,bad.QRS);
    else
        bad.QRS = [];
    end
else
    bad.QRS = [];
end

% Standard deviation in average by marker (evoked potentials)
if isfield(settings,'AVG') & ~isempty(settings.AVG) & isfield(settings,'AVGstd') & ~isempty(settings.AVGstd)
    disp('     detect bad activations by std in average')
    for i = 1:length(settings.AVG.marker)
        tmp = header.events.POS(1,ismember(header.events.TYP,settings.AVG.marker{i})==1);
        if ~isempty(settings.AVG.markerOffset) & settings.AVG.markerOffset > 0
            tmp = tmp - int64(settings.AVG.markerOffset);
        end
        tmp = tmp(1,tmp > 0);
        markerevents{i} = tmp(1,tmp < (size(data,2) - settings.AVG.markerlength)); %#ok<AGROW>
        clearvars tmp
    end
    if isfield(settings.AVG,'combinemarker') & settings.AVG.combinemarker == true
        tmp = [];
        for i = 1:length(markerevents)
            tmp = union(tmp,markerevents{i});
        end
        if size(tmp,1) > 1
            tmp = tmp';
        end
        clearvars markerevents
        markerevents = tmp;
        clearvars tmp
    else
        disp('      only first selected marker is used')
        markerevents = markerevents{1};
    end
    if ~isfield(header,'numauxchannels')
        header.numauxchannels = 0;
    end
    average = lab_average_epochs(data,header,markerevents,settings.AVG,settings.AVG.marker);
    for ch=1:header.numdatachannels
        bad.AVGstd(1,ch) = std(average(ch,:)) / std(average(:));
    end
    bad.average = average;
    if isfield(settings,'AVGmode') & strcmp(settings.AVGmode,'detect bad')
        bad.AVG = find(bad.AVGstd >= settings.AVGstd);
    else
        bad.AVG = find(bad.AVGstd < settings.AVGstd);
    end
    badelectrodes = union(badelectrodes,bad.AVG);
end

% Look for channels without data
eeg_chans = 1:header.numdatachannels;
tmp = zeros(1,header.numdatachannels);
for ch=1:header.numdatachannels
    tmp(ch)=std(data(ch,:));
end
chanNoSignal = tmp<10^-14;
chanNoSignal = eeg_chans(logical(chanNoSignal));
if min(chanNoSignal) > 0
    badelectrodes = union(badelectrodes,chanNoSignal);
    bad.nosignal = chanNoSignal;
else
    bad.nosignal = [];
end
clearvars tmp eeg_chans;

% Looking for jump artifacts
% data = medfilt1(data, 9, [], 2);
% data = abs([diff(data, 1, 2) zeros(size(data,1),1)]);

if isfield(settings,'ecg_chan')
    settings = rmfield(settings,'ecg_chan');
end
if isfield(settings,'blinkchans')
    settings = rmfield(settings,'blinkchans');
end

if size(badelectrodes,1) > 1
    badelectrodes = badelectrodes';
end
bad.percentbad = percentbad;
bad.error = 0;

if isfield(cfg,'EEG_file') & ~exist('noverbose','var')
    [~,~,~,FilenameS_cfg] = lab_filename(cfg.EEG_file);
    Verbose_file=[FilenameS_cfg '_bad.vrb'];
    fid=fopen(fullfile(cfg.EEG_filepath,Verbose_file),'w');
    fprintf(fid,'Detection of bad channels\n');
    fprintf(fid,datestr(now,0));
    fprintf(fid,'\n');
    fprintf(fid,['Input file: ' cfg.EEG_file]);
    fprintf(fid,'\n');
    fprintf(fid,['File-length: ' num2str(header.numtimeframes)]);
    fprintf(fid,'\n\n');
    lab_write_bad_vrb(fid,settings,bad)
    fclose(fid);
end

global dodebug
if isfield(dodebug,'bad_container') & dodebug.bad_container == true
    [~,~,~,FilenameS_mat] = lab_filename(cfg.EEG_file);
    if isfield(dodebug,'bad_container_folder') & exist(dodebug.bad_container_folder,'dir')
        Filepath_mat = dodebug.bad_container_folder;
    else
        Filepath_mat = cfg.EEG_filepath;
    end
    save(fullfile(Filepath_mat,[FilenameS_mat '_bad.mat']),'badelectrodes','bad','settings','-v7.3');
end

return