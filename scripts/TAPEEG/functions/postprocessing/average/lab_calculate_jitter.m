% Function to calculate jitter in averaged erp-signals
%
% Jitter = lab_calculate_jitter(average,epochs,settings)
%
% Written by F. Hatz 2014 Neurology Basel

function Jitter = lab_calculate_jitter(average,epochs,epochsValid,settings)

fprintf('   Calculate Jitter: ');

if ~exist('settings','var') | ~isfield(settings,'range')
    settings = lab_get_JITTER;
    if isempty(settings)
        Jitter = [];
        return
    end
end
if ~exist('epochsValid','var')
    epochsValid = true(size(epochs,1),size(epochs,3));
end

[maxvalue,maxindex] = max(abs(average),[],2);
[~,tmp] = max(maxvalue);
maxvalue = maxvalue(tmp);
tmp1 = mean(abs(average(tmp,1:settings.prerange)));
tmp2 = std(abs(average(tmp,1:settings.prerange)));
SDvalid = ((maxvalue - tmp1) / tmp2) * (settings.percentmax / 100);
for i = 1:size(average,1)
    tmp1 = mean(abs(average(i,1:settings.prerange)));
    tmp2 = std(abs(average(i,1:settings.prerange)));
    tmp = find(abs(average(i,:)) > (tmp1 + SDvalid*tmp2),1,'first');
    if isempty(tmp)
        Jitter.ChanValid(i) = false;
    else
        Jitter.ChanValid(i) = true;
    end
    clearvars tmp tmp1 tmp2
end
chans = 1:size(average,1);
Range = floor(settings.range/2);
Rangej = floor(settings.rangej);

epochsstart = zeros(length(chans),1);
template = zeros(length(chans),2*Range);
epochstmp = zeros(length(chans),2*(Range+Rangej)+1,size(epochs,3));
for i = 1:length(chans)
    if maxindex(i) > Range & (maxindex(i)+Range) < size(average,2)
        template(i,:) = average(chans(i),(maxindex(i)-Range+1):(maxindex(i)+Range));
    else
        template(i,:) = zeros(1,2*Range);
    end
    if (maxindex(i) - (Range+Rangej)) > 0 & (maxindex(i) + Range + Rangej) <= size(epochs,2)
        epochsstart(i,1) = (maxindex(i) - (Range+Rangej));
    elseif (maxindex(i) - (Range+Rangej)) <= 0
        epochsstart(i,1) = 1;
    else
        epochsstart(i,1) = size(epochs,2)-(2*(Range+Rangej)+1)+1;
    end
    epochstmp(i,:,:) = epochs(chans(i),epochsstart(i,1):epochsstart(i,1)+2*(Range+Rangej),:);
end
clearvars i maxindex

Sindex = (size(epochstmp,2) - size(template,2) +1);
Serror = zeros(size(epochstmp,1),(2*Rangej+1),size(epochstmp,3));
counter = round(size(epochs,3) / 20);
template = template - repmat(mean(template,2),[1 size(template,2)]);
for j = 1:size(epochstmp,3)
    if mod(j,counter) == 0
        fprintf('.');
    end
    for i = 1:Sindex
        tmp =epochstmp(:,i:(i+size(template,2)-1),j);
        tmp = tmp - repmat(mean(tmp,2),[1 size(tmp,2)]);
        Serror(:,i,j) = sum((tmp-template).^2,2).^2;
    end
end
disp(':')
clearvars i j counter Sindex tmp

minerror = zeros(size(Serror,1),size(Serror,3));
minerrorI = zeros(size(Serror,1),size(Serror,3));
for j = 1:size(Serror,3)
    for i = 1:size(Serror,1)
        [minerror(i,j),minerrorI(i,j)] = min(Serror(i,:,j));
    end
end
clearvars j Serror


if isfield(settings,'editpre') & settings.editpre == false
    epochsstart = epochsstart - settings.prerange;
end
Jitter.LagAll = minerrorI + repmat(epochsstart,1,size(minerrorI,2)) + Range - 1;
Jitter.Valid = true(size(Jitter.LagAll));
Jitter.Valid(epochsValid == 0) = false;
Jitter.Valid(Jitter.ChanValid == false,:) = false;
Jitter.Valid2 = Jitter.Valid;
Jitter.Lag = Jitter.LagAll;
for i = 1:size(minerror,1)
    Jitter.LagAll(i,minerror(i,:) > (mean(minerror(i,:)) + settings.SD*std(minerror(i,:)))) = epochsstart(i,1) + Range + 2*Rangej;
    Jitter.Valid(i,minerror(i,:) > (mean(minerror(i,:)) + settings.SD*std(minerror(i,:)))) = false;
end
Jitter.Lag(Jitter.Valid == false) = NaN;
Jitter.Mean = nanmean(Jitter.Lag,2);
Jitter.Std = nanstd(Jitter.Lag,[],2);
Jitter.ErrorAll = minerror;
Jitter.Error = Jitter.ErrorAll;
Jitter.Error(Jitter.Valid == false) = NaN;
Jitter.MeanError = nanmean(Jitter.Error,2);

clearvars j minerror minerrorI