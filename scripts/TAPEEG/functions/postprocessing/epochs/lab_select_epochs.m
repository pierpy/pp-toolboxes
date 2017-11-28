% Helper script for lab_save_epochs / Function to select epochs
%
% [epochs,settings] = lab_select_epochs(data,header,settings)
%
% data = matrix (channels x timeframes)
% header = (header.samplingrate is used)
% settings.length = length of epochs
% settings.BAD = settings for bad channels detection
%
% Written by F. Hatz 2012

function [epochs,settings,header] = lab_select_epochs(data,header,settings,length_tot)
disp('   Find valid epochs...')
epochs = [];

if ~exist('settings','var') | ~isfield(settings,'length')
    settings.length = 4096;
end
if ~exist('length_tot','var') | isempty(length_tot) | length_tot < settings.length
    length_tot = settings.length;
end
if isfield(settings,'percentgood')
    percentgood = settings.percentgood / 100;
else
    percentgood = 1;
end
shiftnr = floor(settings.length / (header.samplingrate / 2));
shift = floor(settings.length / shiftnr);
if settings.length < header.samplingrate
    shiftnr = 1;
    shift = settings.length;
elseif shift < header.samplingrate;
    shiftnr = floor(settings.length / header.samplingrate);
    shift = floor(settings.length / shiftnr);
end
if length_tot < header.samplingrate
    shiftnr2 = 1;
else
    shiftnr2 = floor(length_tot / header.samplingrate);
end
epochsnr = floor(size(data,2) / shift);
numtimeframes = epochsnr * shift;
data = data(:,1:numtimeframes);
startflag = 1;
stopflag = size(data,2);

% Search bad shifts
if isfield(header,'bad') & isfield(header.bad,'epochs') & size(header.bad.epochs,2) == epochsnr
    bad = header.bad;
    electrodesValid = abs(bad.epochs - 1);
else
    disp('   Search for bad channels in epochs')
    if ~isfield(settings,'BAD')
        settings.BAD.freqlim50 = 50;
        settings.BAD.freqlim60 = 50;
        settings.BAD.freqlimlow = 50;
        settings.BAD.freqlimhigh = 50;
        settings.BAD.spectslow = [0 2];
        settings.BAD.spectshigh = [15 50];
        settings.BAD.zvaluevars = 3;
        settings.BAD.zvaluehurst = 3;
    end
    settingsB = settings.BAD;
    settingsB.length = shift / header.samplingrate;
    [~,bad] = lab_detect_bad(data,header,settingsB,[],'noverbose');
    header.bad = bad;
    electrodesValid = abs(bad.epochs - 1);
    clearvars settingsB
end

electrodesMarker = ones(size(electrodesValid));
% search bad shifts by markers
if isfield(settings,'markerexclude') & ~isempty(settings.markerexclude) & isfield(header,'events')
    startE = [];
    stopE = [];
    for i = 1:size(settings.markerexclude,1)
        tmp = find(strcmp(header.events.TYP,settings.markerexclude(i,1)));
        if ~isempty(tmp)
            startE = [startE header.events.POS(tmp)]; %#ok<AGROW>
            stopE = [stopE (header.events.POS(tmp)+header.events.DUR(tmp))]; %#ok<AGROW>
        end
    end
    startE = ceil(double(startE) / shift);
    startE(startE > size(electrodesValid,2)) = size(electrodesValid,2);
    startE(startE < 1) = 1;
    stopE = ceil(double(stopE) / shift);
    stopE(stopE > size(electrodesValid,2)) = size(electrodesValid,2);
    stopE(stopE < 1) = 1;
    for i = 1:size(startE,2)
        electrodesValid(:,startE(i):stopE(i)) = 0;
        electrodesMarker(:,startE(i):stopE(i)) = 0;
    end
    clearvars startE stopE
end

% Look for marker start and label all epochs before that marker as bad
if isfield(settings,'markerstart') & ~isempty(settings.markerstart) & isfield(header,'events')
    tmp = find(strcmp(header.events.TYP,settings.markerstart), 1, 'last' );
    tmp = header.events.POS(tmp);
    tmp = ceil(tmp/shift);
    if tmp > size(electrodesValid,2)
        tmp = size(electrodesValid,2);
    elseif tmp < 1
        tmp = 1;
    end
    if ~isempty(tmp)
        electrodesValid(:,1:tmp) = 0;
        electrodesMarker(:,1:tmp) = 0;
        startflag = tmp;
    end
end

% Look for marker stop and label all epochs after that marker as bad
if isfield(settings,'markerstop') & ~isempty(settings.markerstop) & isfield(header,'events')
    tmp = find(strcmp(header.events.TYP,settings.markerstop), 1 );
    tmp = header.events.POS(tmp);
    tmp = ceil(tmp/shift);
    if tmp > size(electrodesValid,2)
        tmp = size(electrodesValid,2);
    elseif tmp < 1
        tmp = 1;
    end
    if ~isempty(tmp)
        electrodesValid(:,tmp:end) = 0;
        electrodesMarker(:,tmp:end) = 0;
        stopflag = tmp;
    end
end

% restrict good epochs to markers
if isfield(settings,'markerinclude') & ~isempty(settings.markerinclude) & ...
        ~isempty(settings.markerinclude{1}) & isfield(header,'events')
    electrodesValid2 = zeros(size(electrodesValid));
    startE = [];
    stopE = [];
    for i = 1:size(settings.markerinclude,1)
        tmp = find(strcmp(header.events.TYP,settings.markerinclude(i,1)));
        if ~isempty(tmp)
            startE = [startE header.events.POS(tmp)]; %#ok<AGROW>
            stopE = [stopE (header.events.POS(tmp)+header.events.DUR(tmp))]; %#ok<AGROW>
        end
        startE = ceil(double(startE) / shift);
        stopE = floor(double(stopE) / shift);
    end
    startE(startE > size(electrodesValid,2)) = size(electrodesValid,2);
    startE(startE < 1) = 1;
    stopE(stopE > size(electrodesValid,2)) = size(electrodesValid,2);
    stopE(stopE < 1) = 1;
    for i = 1:size(startE,2)
        electrodesValid2(:,startE(i):stopE(i)) = 1;
    end
    electrodesValid(~electrodesValid2) = 0;
    clearvars electrodesValid2 startE stopE
end

if isfield(bad,'peak2min') & ~isempty(bad.peak2min)
    electrodesValid = electrodesValid .* repmat(bad.peak2min,size(electrodesValid,1),1);
end
if isfield(bad,'microcorr') & ~isempty(bad.microcorr)
    electrodesValid = electrodesValid .* repmat(bad.microcorr,size(electrodesValid,1),1);
end

% Calculate mean for epochs
for i = 1:(epochsnr - shiftnr2 + 1)
    epochsVmean(:,i) = mean(electrodesValid(:,i:(i+shiftnr2-1)),2); %#ok<AGROW>
    epochsSV.Marker(:,i) = mean(electrodesMarker(:,i:(i+shiftnr2-1)),2);
    if isfield(bad,'epochsSpect50') & ~isempty(bad.epochsSpect50)
        epochsSV.Freq50(:,i) = mean(abs(bad.epochsSpect50(:,i:(i+shiftnr2-1))-1),2);
    end
    if isfield(bad,'epochsSpect60') & ~isempty(bad.epochsSpect60)
        epochsSV.Freq60(:,i) = mean(abs(bad.epochsSpect60(:,i:(i+shiftnr2-1))-1),2);
    end
    if isfield(bad,'epochsSpectLow') & ~isempty(bad.epochsSpectLow)
        epochsSV.FreqLow(:,i) = mean(abs(bad.epochsSpectLow(:,i:(i+shiftnr2-1))-1),2);
    end
    if isfield(bad,'epochsSpectHigh') & ~isempty(bad.epochsSpectHigh)
        epochsSV.FreqHigh(:,i) = mean(abs(bad.epochsSpectHigh(:,i:(i+shiftnr2-1))-1),2);
    end
    if isfield(bad,'epochsVar') & ~isempty(bad.epochsVar)
        epochsSV.Var(:,i) = mean(abs(bad.epochsVar(:,i:(i+shiftnr2-1))-1),2);
    end
    if isfield(bad,'epochsHurst') & ~isempty(bad.epochsHurst)
        epochsSV.Hurst(:,i) = mean(abs(bad.epochsHurst(:,i:(i+shiftnr2-1))-1),2);
    end
    if isfield(bad,'peak2min') & ~isempty(bad.peak2min)
        epochsSV.peak2min(1,i) = mean(bad.peak2min(1,i:(i+shiftnr2-1)),2);
    end
    if isfield(bad,'microcorr') & ~isempty(bad.microcorr)
        epochsSV.microcorr(1,i) = mean(bad.microcorr(1,i:(i+shiftnr2-1)),2);
    end
end

% set minimalpart to integer value = number of epochs
if ~isfield(settings,'minimalpart')
    settings.minimalpart = 0.7;
end
if settings.minimalpart < 1
    minimalpart = round(settings.minimalpart * ((stopflag - startflag + 1) / settings.length));
    if minimalpart < 1
        minimalpart = 1;
    end
elseif round(settings.minimalpart)*settings.length > (stopflag - startflag + 1)
    minimalpart = floor((stopflag - startflag + 1)/settings.length);
else
    minimalpart = round(settings.minimalpart);
end

% find good epochs (sensitivity is increased until minimalpart is reached)
goodEpochs = [];
while size(goodEpochs,2) < minimalpart
    goodEpochs = [];
    while isempty(goodEpochs)
        goodEpochs = mean(epochsVmean,1);
        goodEpochs = find(goodEpochs >= percentgood);
        if isempty(goodEpochs)
            percentgood = percentgood * 0.99;
        end
        if percentgood < 0.01
            epochs = [];
            return
        end
    end
    i = 1;
    while i < size(goodEpochs,2)
        if (goodEpochs(1,i) + shiftnr) >= goodEpochs(1,(i+1))
            if (i+1) < size(goodEpochs,2)
                goodEpochs = [goodEpochs(1,1:i) goodEpochs(1,(i+2):end)];
            end
        end
        if (goodEpochs(1,i) + shiftnr) <= goodEpochs(1,(i+1)) | i == (size(goodEpochs,2) - 1)
            i = i + 1;
        end
    end
    if size(goodEpochs,2) > 1 & (goodEpochs(1,end -1) + shiftnr2) >= goodEpochs(1,end)
        goodEpochs = goodEpochs(1,1:(end-1));
    end
    if size(goodEpochs,2) < minimalpart
        percentgood = percentgood * 0.99;
    end
    if percentgood < 0.01
        break
    end
end

% sort epochs by quality
goodsum = zeros(size(goodEpochs));
badoutelectrodes = zeros(size(goodEpochs,2),size(electrodesValid,1));
for i = 1:size(goodEpochs,2)
    tmp = electrodesValid(:,goodEpochs(1,i):(goodEpochs(1,i)+shiftnr2-1));
    goodsum(i) = sum(tmp(:))/length(tmp(:));
    tmp = 1-tmp;
    badoutelectrodes(i,:) = sum(tmp,2) / shiftnr2;
    clearvars tmp
end
[goodsum,goodsort] = sort(goodsum,'descend');
if length(goodsort) > minimalpart
    goodsort = goodsort(1:minimalpart);
end
goodEpochs = goodEpochs(goodsort);
badoutelectrodes = badoutelectrodes(goodsort,:);

for i = 1:size(goodEpochs,2)
    marker = (goodEpochs(1,i)-1)*shift;
    if (marker + length_tot) > size(data,2)
        marker = size(data,2) - length_tot;
    end
    epochs.data(:,:,i) = detrend(data(:,(marker +1):(marker + length_tot))')';
    epochs.badchans{1,i} = find(epochsVmean(:,goodEpochs(1,i)) < percentgood)';
    if isfield(epochsSV,'Freq50')
        epochs.bad.Spect50{1,i} = find(epochsSV.Freq50(:,goodEpochs(1,i)) < percentgood)';
    end
    if isfield(epochsSV,'Freq60')
        epochs.bad.Spect60{1,i} = find(epochsSV.Freq60(:,goodEpochs(1,i)) < percentgood)';
    end
    if isfield(epochsSV,'FreqLow')
        epochs.bad.SpectLow{1,i} = find(epochsSV.FreqLow(:,goodEpochs(1,i)) < percentgood)';
    end
    if isfield(epochsSV,'FreqHigh')
        epochs.bad.SpectHigh{1,i} = find(epochsSV.FreqHigh(:,goodEpochs(1,i)) < percentgood)';
    end
    if isfield(epochsSV,'Var')
        epochs.bad.Var{1,i} = find(epochsSV.Var(:,goodEpochs(1,i)) < percentgood)';
    end
    if isfield(epochsSV,'Hurst')
        epochs.bad.Hurst{1,i} = find(epochsSV.Hurst(:,goodEpochs(1,i)) < percentgood)';
    end
    if isfield(epochsSV,'Marker')
        epochs.bad.Marker{1,i} = find(epochsSV.Marker(:,goodEpochs(1,i)) < percentgood)';
    end
    if isfield(epochsSV,'peak2min')
        epochs.bad.Peak2Min(1,i) = epochsSV.peak2min(1,goodEpochs(1,i));
    end
    epochs.events.POS(1,i) = int64(marker+1);
    epochs.events.DUR(1,i) = int64(length_tot);
    epochs.events.OFF(1,i) = int64(0);
    epochs.events.TYP{1,i} = ['Epoch' num2str(i)];
end
[~,tmp] = sort(epochs.events.POS);
epochs.events.POS = epochs.events.POS(1,tmp);
epochs.events.DUR = epochs.events.DUR(1,tmp);
epochs.events.OFF = epochs.events.OFF(1,tmp);
epochs.events.TYP = epochs.events.TYP(1,tmp);
clearvars tmp

epochs.goodsum = goodsum;
epochs.badoutelectrodes = badoutelectrodes;
epochs.goodmarkers = goodEpochs; 
epochs.markersvalid = electrodesValid;
epochs.markersvalidshift = shiftnr;
epochs.markersvalidlength = shift;
epochs.settings = settings.BAD;
epochs.settings.minimalpart = minimalpart;
epochs.percentgood = percentgood * 100;
if isfield(settings,'interpolatebad')
    epochs.interpolatebad = settings.interpolatebad;
end

return


