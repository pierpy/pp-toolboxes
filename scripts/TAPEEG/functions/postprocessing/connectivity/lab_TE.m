function [TE,TE_all,MI,MI_all,TE_delays] = lab_TE(data,header,settings,cfg)
    
global TEraw

if ~exist('cfg','var')
    cfg = [];
end
if ~exist('settings','var') | ~isfield(settings,'el')
    settings.el = 5;
    settings.tau = 0.35;
    settings.delay_min = 1;
    settings.delay_step = 2;
    settings.delay_max = 50;
end

if ~isfield(cfg,'firstfilefolder')
    TEraw = {data};
elseif cfg.firstfilefolder
    TEraw = {data};
else
    TEraw{end+1} = data;
end

if isfield(cfg,'lastfilefolder') & cfg.lastfilefolder == false
    TE = [];
    TE_all = {};
    TE_delays = [];
    return
elseif isfield(cfg,'lastsegment') & cfg.lastsegment == false
    TE = [];
    TE_all = {};
    TE_delays = [];
    return
end

clear data
for i = 1:length(TEraw);
    data.trial{i} = TEraw{i};
    data.time{i} = (1:header.numtimeframes)./header.samplingrate; 
end
for i=1:size(data.trial{1},1)
    data.label{i}=num2str(i);
end
data.fsample = header.samplingrate;
[A,B] = meshgrid(1:size(data.trial{1},1),1:size(data.trial{1},1));
c=cat(2,A',B');
d=reshape(c,[],2);d(d(:,1)==d(:,2),:)=[];
clear A B c

data.TEprepare.nrtrials = ones(size(d))*size(data.trial,2);

data.TEprepare.channelcombi = d;
data.TEprepare.channelcombilabel = data.label(d);

TEcfg = [];
TEcfg.sgncmb = data.TEprepare.channelcombi;
TEcfg.channelcombilabel = data.label(d);
clear d

data.TEprepare.trials = {}; 
for i=1:size(TEcfg.sgncmb,1)
    for j=1:size(TEcfg.sgncmb,2)
        data.TEprepare.trials{i,j} = 1:length(TEraw); 
    end
end
data.TEprepare.timeindices = [1 size(data.trial{1},2)];
TEcfg.toi                 = [min(data.time{1}),max(data.time{1})]; % time of interest

% estimator
TEcfg.TEcalctype  = 'VW_ds'; % use the new TE estimator (Wibral, 2013)

% kernel-based TE estimation
%TEcfg.flagNei = 'Mass' ;           % neigbour analyse type
TEcfg.kth_neighbors = 4;   
% calculate ACT
data.TEprepare.ACT = TEactdetect(data2datacell(data.TEprepare.channelcombi,data),1000,data.TEprepare.timeindices);
TEcfg.TheilerT = 'ACT';
TEcfg.calctime = 'no';
TEcfg.trialselect = 'No';
TEcfg.extracond = 'none';

prediction_times = settings.delay_min:settings.delay_step:settings.delay_max;      % minimum u to be scanned

TEcfg.dim = ones(1,size(data.TEprepare.trials,1))*settings.el;
TEcfg.tau = ones(1,size(data.TEprepare.trials,1))*settings.tau;

warning off %#ok<WNOFF>
for i=1:length(prediction_times)
    data.TEprepare.u_in_samples = prediction_times(i);
    TEcfg.predicttime_u = prediction_times(i);
    [TEresult(i)] = transferentropy(TEcfg,data); %#ok<AGROW>
    TEtmp = mean(TEresult(i).TEmat,2);
    TEtmp = reshape(TEtmp,max(TEresult(i).sgncmb(:))-1,max(TEresult(i).sgncmb(:)));
    TE_all(:,:,i) = lab_add_diagonal(TEtmp); %#ok<AGROW>
    MItmp = mean(TEresult(i).MImat,2);
    MItmp = reshape(MItmp,max(TEresult(i).sgncmb(:))-1,max(TEresult(i).sgncmb(:)));
    MI_all(:,:,i) = lab_add_diagonal(MItmp); %#ok<AGROW>
end
warning on %#ok<WNON>
TE = max(TE_all,[],3);
MI = max(MI_all,[],3);
TE_delays = prediction_times;

return

function datacell = data2datacell(channelcombi,data)
% create datacell {channelcombi x 2} including the matrix (trial x
% timepoints) for each channel.
datacell = cell(size(channelcombi,1),2);
for cc = 1:size(channelcombi,1)
    for pp = 1:2
        datamat = zeros(size(data.trial,2),size(data.trial{1},2));
        for ii = 1:size(data.trial,2)
            datamat(ii,:)=data.trial{ii}(channelcombi(cc,pp),:);
        end
        datacell{cc,pp}=datamat;
        clear datamat;
    end
end