function Phase = lab_calculate_connect_phaseepochs(phase,header,settings,Epochs,Idx)

Phase = [];
shift = Epochs.markersvalidlength;
shiftnr = Epochs.markersvalidshift;
electrodesValid = Epochs.markersvalid;
length_tot = size(Epochs.data,2);
percentgood = Epochs.percentgood / 100;
if length_tot < header.samplingrate
    shiftnr2 = 1;
else
    shiftnr2 = floor(length_tot / header.samplingrate);
end

electrodesValid = repmat(electrodesValid,[1 1 shift]);
electrodesValid = permute(electrodesValid,[1 3 2]);
electrodesValid = reshape(electrodesValid,[size(electrodesValid,1) size(electrodesValid,2)*size(electrodesValid,3)]);
electrodesValid = electrodesValid(:,Idx(Idx<=size(electrodesValid,2)));

epochsnr = floor(size(electrodesValid,2) / shift);
electrodesValid = electrodesValid(:,1:epochsnr*shift);
epochsVmean = reshape(electrodesValid,[size(electrodesValid,1) shift epochsnr]);
epochsVmean = mean(epochsVmean,2);
epochsVmean = permute(epochsVmean,[1 3 2]);

% set minimalpart to integer value = number of epochs
if ~isfield(settings,'minimalpart')
    settings.minimalpart = 0.7;
end
if settings.minimalpart < 1
    minimalpart = round(settings.minimalpart * (size(electrodesValid,2) / settings.length));
    if minimalpart < 1
        minimalpart = 1;
    end
elseif round(settings.minimalpart)*settings.length > size(electrodesValid,2)
    minimalpart = floor(size(electrodesValid,2)/settings.length);
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
            Phase = [];
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
[~,goodsort] = sort(goodsum,'descend');
if length(goodsort) > minimalpart
    goodsort = goodsort(1:minimalpart);
end
goodEpochs = goodEpochs(goodsort);

for i = 1:size(goodEpochs,2)
    marker = (goodEpochs(1,i)-1)*shift;
    if (marker + length_tot) > size(phase,2)
        marker = size(phase,2) - length_tot;
    end
    Phase(:,:,i) = detrend(phase(:,(marker +1):(marker + length_tot))')'; %#ok<AGROW>
end

return