function res = RaguMSStats(d,h)

% Copyright 2009-2011 Thomas Koenig
% distributed under the terms of the GNU AFFERO General Public License

nClas = size(d.MSMaps,1);
%rng('shuffle');
%s = RandStream.create('mt19937ar','seed',sum(100*clock));

%if verLessThan('matlab', '7.0.1')
%    rand('seed',sum(100*clock));
%else
%    RandStream.setGlobalStream(s);
%end

SelIndFeature = d.IndFeature;
SelIndFeature(isnan(d.IndFeature)) = [];
SelDesign = d.Design;
SelDesign(isnan(SelDesign(:,1)),:) = [];

DoGroup = (numel(unique(SelIndFeature))> 1);
DoF1    = (numel(unique(SelDesign(:,1)))> 1);

if d.ContF1 == 1
    errordlg('Continuous predictors currently not supported');
end

[NewDesign,iDesign,CondIndices,DoF2] = Ragu_SortOutWithinDesign(d.Design',d.TwoFactors);

nFeatures = 6;
d.MSEffectSize  = zeros(2,4,nClas,d.Iterations,nFeatures);
if ~isfield(d,'MDDelta')
    d.MDDelta       = zeros(2,4,d.Iterations,15);
end

if nargin < 2
    h = waitbar(0,'Computing microstate stats, please wait...');
end

set(h,'WindowStyle','modal');


InData = zeros(size(d.V,1),numel(iDesign),size(d.V,3),size(d.V,4));

for i = 1:numel(iDesign)
    InData(:,i,:,:) = mean(d.V(:,CondIndices == i,:,:),2);
end

if (DoF1 && DoF2)
    if abs(corr(NewDesign(:,1),NewDesign(:,2))) > 0.01
        error('Design not orthogonal');
    end
end

%InData = d.V;

%InData(:,isnan(d.Design(:,1)),:,:) = [];
InData(isnan(d.IndFeature),:,:,:) = [];

tic

[EffectSize,dlt] = AllMSEffectSizes(InData,d.MSMaps,NewDesign',iDesign,SelIndFeature,DoGroup,DoF1,DoF2,d.bSmoothLabels,d.nWindowSize,d.LabelPenalty,d.StartFrame,d.EndFrame);

d.MSEffectSize(:,:,:,1,:) = EffectSize;

d.MDDelta(:,:,1,nClas) = dlt;

%if nargin < 2
    waitbar(1/d.Iterations,h);
    set(h,'Name',sprintf('Remaining time: %01.0f:%02.0f min',floor(toc()*(d.Iterations)/60),rem(toc()*(d.Iterations-1),60)));

%else            
%ShowProgress(1/d.Iterations,h);
%end

NoXing = d.NoXing;
Iterations = d.Iterations;

for iter = 2:d.Iterations
    r_in = zeros(size(InData));    
    for j = 1:size(InData,1)
        r_in(j,:,:) = InData(j,PermDesign(NewDesign',NoXing),:);
    end
    rIndFeat = SelIndFeature(randperm(numel(SelIndFeature)));

%    [On,Off,Dur,auc,cog,m_gfp,dlt] = AllMSEffectSizes(r_in,d.MSMaps,uDesign',iDesign,rIndFeat,DoGroup,DoF1,DoF2,d.bSmoothLabels,d.nWindowSize,d.LabelPenalty,d.StartFrame,d.EndFrame);
    [EffectSize,dlt] = AllMSEffectSizes(r_in,d.MSMaps,NewDesign',iDesign,rIndFeat,DoGroup,DoF1,DoF2,d.bSmoothLabels,d.nWindowSize,d.LabelPenalty,d.StartFrame,d.EndFrame);

    d.MSEffectSize(:,:,:,iter,:) = EffectSize;


    d.MDDelta(:,:,iter,nClas) = dlt;
    if nargin < 2
        waitbar(iter/Iterations,h);
        set(h,'Name',sprintf('Remaining time: %01.0f:%02.0f min',floor(toc()*(Iterations/iter-1)/60),rem(toc()*(Iterations/iter-1),60)));
%    else
%         ShowProgress(iter/Iterations,h);
    end
end

d.pMSStats = nan(size(d.MSEffectSize));

for i = 1:2
    for j = 1:4
        for c = 1:nClas
            for f = 1:6
                if max(squeeze(d.MSEffectSize(i,j,c,:,f))) == 0
                    d.pMSStats(i,j,c,:,f) = 1;
                else
                    inp = squeeze(d.MSEffectSize(i,j,c,:,f));
                    idx = find(~isnan(inp));
                    for r = 1:numel(inp)
                        if isnan(inp(r))
                            d.pMSStats(i,j,c,r,f) = NaN;
                        else
                            d.pMSStats(i,j,c,r,f) = sum(inp(r) <= inp(idx)) / numel(idx);
                        end
                    end
                end
            end
        end
    end
end

OverallDelta = squeeze(sum(sum(d.MDDelta(:,:,:,nClas),2),1));
Excess = OverallDelta(1) - mean(OverallDelta(2:end));
OverallP = sum(OverallDelta(1) <= OverallDelta) / numel(OverallDelta);

if (nargin < 2)
    close(h);
end
res = d;

%function [On,Off,Dur,auc,cog,gfp,dlt] = AllMSEffectSizes(Data,MSMaps,Design,iDesign,IndFeat,DoGroup,DoF1,DoF2,bSmooth,b,p,StartFrame,EndFrame)
function [EffectSize,dlt] = AllMSEffectSizes(Data,MSMaps,Design,iDesign,IndFeat,DoGroup,DoF1,DoF2,bSmooth,b,p,StartFrame,EndFrame)

nClas = size(MSMaps,1);
nFeatures = 6;
EffectSize = zeros(2,4,nClas,nFeatures);

dlt = zeros(2,4);

AllC = ones(size(Design,1),1);
AllC(isnan(Design(:,1)),1) = nan;
AllG = ones(size(IndFeat));
AllG(isnan(IndFeat)) = nan;

%Group main effect
if DoGroup
    gm = RaguGrandMeans(Data,IndFeat,AllC);
    [MSClass,MSFit,MSgfp] = RaguFitMicrostates(gm,MSMaps,bSmooth,b,p,StartFrame,EndFrame);

    [GMEffect,ndtGM] = RaguMSFeatures(MSClass,nClas,MSFit,MSgfp);
    
    for i = 1:nFeatures
        EffectSize(2,1,:,i) = MSvar(GMEffect(:,:,:,i));
    end

    dlt(2,1,:) = ndtGM;
else
    GMEffect = zeros(1,1,1,nFeatures);

    ndtGM = 0;
end
    
% F1  
if DoF1
    gm = RaguGrandMeans(Data,AllG,Design(:,1));
    [MSClass,MSFit,MSgfp] = RaguFitMicrostates(gm,MSMaps,bSmooth,b,p,StartFrame,EndFrame);
    
    [F1Effect,ndtF1] = RaguMSFeatures(MSClass,nClas,MSFit,MSgfp);

    for i = 1:nFeatures
        EffectSize(1,2,:,i) = MSvar(F1Effect(:,:,:,i));
    end
    dlt(1,2)   = ndtF1;
else
    F1Effect = zeros(1,1,1,nFeatures);

    ndtF1 = 0;
end

%F1 x G
if DoGroup && DoF1
    gm = RaguGrandMeans(Data,IndFeat,Design(:,1));
    [MSClass,MSFit,MSgfp] = RaguFitMicrostates(gm,MSMaps,bSmooth,b,p,StartFrame,EndFrame);

    [F1GEffect,ndtF1G] = RaguMSFeatures(MSClass,nClas,MSFit,MSgfp);    
    for i = 1:nFeatures
        EffectSize(2,2,:,i) = MSvar(F1GEffect(:,:,:,i)  - ExpandResults(GMEffect(:,:,:,i) ,AllC,IndFeat,Design(:,1),IndFeat) - ExpandResults(F1Effect(:,:,:,i) ,Design(:,1),AllG,Design(:,1),IndFeat));
    end
    dlt(2,2  ) = ndtF1G; % - ndtGM - ndtF1;
else
    F1GEffect = zeros(1,1,1,nFeatures);

    ndtF1G = 0;
end

% F2
if DoF2
    gm = RaguGrandMeans(Data,AllG,Design(:,2));
    [MSClass,MSFit,MSgfp] = RaguFitMicrostates(gm,MSMaps,bSmooth,b,p,StartFrame,EndFrame);

    [F2Effect,ndtF2] = RaguMSFeatures(MSClass,nClas,MSFit,MSgfp);

    for i = 1:nFeatures
        EffectSize(1,3,:,i) = MSvar(F2Effect(:,:,:,i));
    end
    
    dlt(1,3  ) = ndtF2;
else
    F2Effect = zeros(1,nFeatures);

    ndtF2 = 0;
end
    
%F2 x G
if DoGroup && DoF2
    gm = RaguGrandMeans(Data,IndFeat,Design(:,2));
    [MSClass,MSFit,MSgfp] = RaguFitMicrostates(gm,MSMaps,bSmooth,b,p,StartFrame,EndFrame);
    
    [F2GEffect,ndtF2G] = RaguMSFeatures(MSClass,nClas,MSFit,MSgfp);
    
    for i = 1:nFeatures
        EffectSize(2,3,:,i) = MSvar(F2GEffect(:,:,:,i)  - ExpandResults(GMEffect(:,:,:,i) ,AllC,IndFeat,Design(:,2),IndFeat) - ExpandResults(F2Effect(:,:,:,i) ,Design(:,1),AllG,Design(:,2),IndFeat));
    end
    
else
    F2GEffect = zeros(1,1,1,nFeatures);
    ndtF2G = 0;
end
    
    
% F1 x F2
if DoF1 && DoF2
%    gm = RaguGrandMeans(Data,AllG,iDesign);
    gm = RaguGrandMeans(Data,AllG,Design);  % TK 3.Mai 2011
    [MSClass,MSFit,MSgfp] = RaguFitMicrostates(gm,MSMaps,bSmooth,b,p,StartFrame,EndFrame);

    [F1F2Effect,ndtF1F2] = RaguMSFeatures(MSClass,nClas,MSFit,MSgfp);

    for i = 1:nFeatures
        EffectSize(1,4,:,i) = MSvar(F1F2Effect(:,:,:,i)  - ExpandResults(F1Effect(:,:,:,i) ,Design(:,1),AllG,iDesign,AllG)  - ExpandResults(F2Effect(:,:,:,i) ,Design(:,2),AllG,iDesign,AllG));
    end
    dlt(1,4  ) = ndtF1F2; % - ndtF1 - ndtF2;
else
    F1F2Effect = zeros(1,1,1,nFeatures);

    ndtF1F2 = 0;
end

% F1 x F2 * G
if DoF1 && DoF2 && DoGroup
%    gm = RaguGrandMeans(Data,IndFeat,iDesign);
    gm = RaguGrandMeans(Data,IndFeat,Design);   % TK 3. Mai 2011
    [MSClass,MSFit,MSgfp] = RaguFitMicrostates(gm,MSMaps,bSmooth,b,p,StartFrame,EndFrame);

    [F1F2GEffect,ndtF1F2G] = RaguMSFeatures(MSClass,nClas,MSFit,MSgfp);

    for i = 1:nFeatures
        EffectSize(2,4,:,i)  = MSvar(F1F2GEffect(:,:,:,i)  - ExpandResults(F1Effect(:,:,:,i) ,Design(:,1),AllG,iDesign,IndFeat)  - ExpandResults(F2Effect(:,:,:,i) ,Design(:,2),AllG,iDesign,IndFeat) - ExpandResults(F1GEffect(:,:,:,i) ,Design(:,1),IndFeat,iDesign,IndFeat)  - ExpandResults(F2GEffect(:,:,:,i) ,Design(:,2),IndFeat,iDesign,IndFeat) - ExpandResults(F1F2Effect(:,:,:,i) ,iDesign,AllG,iDesign,IndFeat));
    end
    dlt(2,4  ) = ndtF1F2G; % - ndtF1 - ndtF2 - ndtF1G - ndtF2G - ndtF1F2;
end

function res = MSvar(in)

    inr = reshape(in,size(in,1) * size(in,2),size(in,3));
    res = var(inr,0,1);


function res = ExpandResults(in,InDesign,InGroup,OutDesign,OutGroup)

IndResultMatrix = zeros(numel(InGroup),size(InDesign,1),size(in,3));

[cii,cji,cki] = unique(InDesign);
[gii,gji,gki] = unique(InGroup);

for c = 1:numel(cii)
    c_idx = (cki == c);
    for g = 1:numel(gii);
        g_idx = (gki == g);
        for n = 1:size(in,3)
            IndResultMatrix(g_idx,c_idx,n) = in(g,c,n);
        end
    end
end

[cii,cji,cki] = unique(OutDesign);
[gii,gji,gki] = unique(OutGroup);

res = zeros(numel(gii),numel(cii),size(in,3));

for c = 1:numel(cii)
    c_idx = (cki == c);
    for g = 1:numel(gii);
        g_idx = (gki == g);
        for n = 1:size(in,3)
            val = IndResultMatrix(g_idx,c_idx,n);
            res(g,c,n) = mean(val(:));
        end
    end
end