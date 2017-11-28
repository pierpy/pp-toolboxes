% Function to average data by markers
%
% [average,epochs,epochsValid] = lab_average_epochs(data,header,markers,settings)
%
% data = matrix (channels x timeframes)
% header = see 'lab_create header'
% markers = vector (1 x eventpositions)
% settings = see lab_set_average
%
% Written by F. Hatz 2015 Neurology Basel

function [average,epochs,epochsValid,epochsBad] = lab_average_epochs(data,header,markers,settings,markerexclude)

if ~exist('markerexclude','var')
    markerexclude = {};
end
if ~exist('settings','var') | isempty(settings) | ~exist('markers','var') | isempty(markers)
    average = [];
    epochs = [];
    epochsValid = [];
    return
end
if isfield(settings,'markerOffset') & settings.markerOffset > 0 & ...
        isfield(settings,'correctbaseline') & strcmp(settings.correctbaseline,'pre-marker period')
    fprintf(['   Averaging including baseline-correction using first ' num2str(settings.markerOffset) ' TFs:'])
    correctbaseline = settings.markerOffset;
elseif isfield(settings,'correctbaseline') & strcmp(settings.correctbaseline,'all TFs')
    fprintf('   Averaging including baseline-correction using all TFs:')
    correctbaseline = settings.markerlength;
else
    fprintf('   Averaging without baseline-correction:')
    correctbaseline = [];
end

epochs = zeros(size(data,1),settings.markerlength,length(markers));
Ptime = 1:size(data,2);
timeline = [];
NumE = length(markers);
CounterE = floor(NumE / 20);
for i = 1:NumE
    if ~isempty(correctbaseline)
        tmp = data(:,markers(i):(markers(i)+settings.markerlength-1));
        epochs(:,:,i) = tmp - repmat(mean(tmp(:,1:correctbaseline),2),1,size(tmp,2));
    else
        epochs(:,:,i) = data(:,markers(i):(markers(i)+settings.markerlength-1));
    end
    timeline = cat(2,timeline,Ptime(1,markers(i):(markers(i)+settings.markerlength-1)));
    if mod(i,CounterE) == 0
        fprintf('.')
    end
end
disp(':')
headerepochs = lab_reallocate_markers(header,timeline,markerexclude);
clearvars tmp

% Reject epochs (FASTER or TAPEEG routines)
epochsBad = [];
if isfield(settings,'Reject') & isfield(settings.Reject,'routine') & strcmp(settings.Reject.routine,'TAPEEG') & ...
        isstruct(settings.Reject.BAD) & ~isempty(settings.Reject.BAD)
    settings.Reject.BAD.length = settings.markerlength / header.samplingrate;
    [~,epochsBad,settings.Reject.BAD] = lab_detect_bad(reshape(epochs,size(epochs,1),size(epochs,2)*size(epochs,3)),headerepochs,settings.Reject.BAD);
    epochsValid = 1 - epochsBad.epochs;
elseif isfield(settings,'Reject') & isfield(settings.Reject,'routine') & strcmp(settings.Reject.routine,'FASTER') & ...
        isnumeric(settings.Reject.BAD) & ~isempty(settings.Reject.BAD)
    epochsValid = ones(size(data,1),length(markers));
    Counter = floor(size(data,1) / 10);
    for j = 1:size(epochs,1)
        xx = reshape(epochs(j,:,:),1,size(epochs,2)*size(epochs,3));
        maxSD = settings.Reject.BAD * std(xx);
        eegvar = zeros(1,length(markers));
        eegamp = zeros(1,length(markers));
        for i = 1:length(markers)
            eegvar(i) = var(xx(:,((i-1)*settings.markerlength + 1:i*settings.markerlength)));
            eegamp(i) = max(xx(:,((i-1)*settings.markerlength + 1:i*settings.markerlength))) - ...
                min(xx(:,((i-1)*settings.markerlength + 1:i*settings.markerlength)));
        end
        eegamp = mean(eegamp) + settings.Reject.BAD*std(eegamp);
        eegvar = mean(eegvar) + settings.Reject.BAD*std(eegvar);
        for i = 1:length(markers)
            if isfield(header,'numauxchannels') & ~isempty(header.numauxchannels) & j > (size(data,1) - header.numauxchannels)
                epochsValid(j,i) = 1;
            elseif ((std(epochs(j,:,i)) >= maxSD) | (var(epochs(j,:,i)) >= eegvar) | ...
                    ((max(epochs(j,:,i)) - min(epochs(j,:,i))) >= eegamp));
                epochsValid(j,i) = 0;
            end
        end
        if mod(j,Counter) == 0
            fprintf('.');
        end
    end
    disp(':')
else
    epochsValid = ones(size(data,1),length(markers));
end

% if rejection method is set to percent channels, modify epochsValid accordingly
if isfield(settings,'Reject') & isfield(settings.Reject,'method') & strcmp(settings.Reject.method,'percent') & ...
        isfield(settings.Reject,'percent') & ~isempty(settings.Reject.percent)
    tmp = sum(epochsValid,1) ./ size(epochsValid,1);
    idx = zeros(1,length(tmp));
    idx(tmp >= (1-settings.Reject.percent/100)) = 1;
    epochsValid = zeros(size(epochsValid));
    epochsValid(:,idx==1) = 1;
    clearvars tmp idx
end

% Calculate average
average = zeros(size(epochs,1),size(epochs,2));
fprintf(['   Averaging using ' settings.AVGmethod])
Counter = floor(size(data,1) / 10);
for i = 1:size(data,1)
    if sum(epochsValid(i,:) == 1) > 1
        if strcmp(settings.AVGmethod,'median')
            average(i,:) = median(epochs(i,:,epochsValid(i,:) == 1),3);
        else
            average(i,:) = mean(epochs(i,:,epochsValid(i,:) == 1),3);
        end
    elseif sum(epochsValid(i,:) == 1) == 1
        average(i,:) = epochs(i,:,epochsValid(i,:) == 1);
    end
    if mod(i,Counter) == 0
        fprintf('.');
    end
end
disp(':')
