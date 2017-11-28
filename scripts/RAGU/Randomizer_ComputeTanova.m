function result = Randomizer_ComputeTanova(d,h)

% Copyright 2009-2011 Thomas Koenig
% distributed under the terms of the GNU AFFERO General Public License

switch d.Normalize
    case 1
        PreNormalizeMode = false;
        bNormalize = false;
    case 2
        PreNormalizeMode = true;
        bNormalize = false;
    case 3
        errordlg('Dissimilarity mode currently not available','Ragu');
        result = NaN;
        return;
        
%        PreNormalizeMode = false;
%        bNormalize = true;
end
SelIndFeature = d.IndFeature;
SelIndFeature(isnan(d.IndFeature)) = [];
SelDesign = d.Design;

SelDesign(isnan(SelDesign(:,1)),:) = [];

DoGroup = (numel(unique(SelIndFeature))> 1);
DoF1    = (numel(unique(SelDesign(:,1)))> 1);
%DoF2    = (numel(unique(SelDesign(:,2)))> 1);

if d.ContF1 == 1
    DoF2 = 0;
    NewDesign = SelDesign(:,1)';
    CondIndices = find(~isnan(NewDesign));
    iDesign = NewDesign;
else
    [NewDesign,iDesign,CondIndices,DoF2] = Ragu_SortOutWithinDesign(d.Design',d.TwoFactors);
end

if d.MeanInterval
    TanovaEffectSize = ones(2,4,1,d.Iterations);
else
    TanovaEffectSize = ones(2,4,d.EndFrame-d.StartFrame+1,d.Iterations);
end

if nargin < 2
    h = waitbar(0,'Computing, please wait...');
end


if (DoF1 && DoF2)
    if abs(corrcoef(NewDesign(:,1),NewDesign(:,2))) > 0.01
        error('Design not orthogonal');
    end
end

InData = zeros(size(d.V,1),numel(iDesign),size(d.V,3),size(d.V,4));

for i = 1:numel(iDesign)
    InData(:,i,:,:) = mean(d.V(:,CondIndices == i,:,:),2);
end

if d.DoGFP    
    InData = sqrt(mean(InData.^2,3));
end
    
%InData(:,isnan(d.Design(:,1)),:,:) = [];
InData(isnan(d.IndFeature),:,:,:) = [];

if d.MeanInterval
    InData = mean(InData(:,:,:,d.StartFrame:d.EndFrame),4);
else
    InData = InData(:,:,:,d.StartFrame:d.EndFrame);
end

if ~d.DoGFP && PreNormalizeMode
    InData = NormDimL2(InData,3);
end

tic
if d.MeanInterval
    TanovaEffectSize(:,:,1,1) = DoAllEffectSizes(InData,NewDesign',iDesign,SelIndFeature,DoGroup,DoF1,DoF2,bNormalize,d.ContBetween,d.ContF1);
else
    TanovaEffectSize(:,:,:,1) = DoAllEffectSizes(InData,NewDesign',iDesign,SelIndFeature,DoGroup,DoF1,DoF2,bNormalize,d.ContBetween,d.ContF1);
end

if nargin > 1
    close(h);
end

% parfor ii=1:N
%     pause( 0.1 );
%     ppm.increment();
% end


NoXing = d.NoXing;
MeanInterval = d.MeanInterval;
% Iterations = d.Iterations;

if 0
    ppm = ParforProgMon( 'Progress',d.Iterations);
    parfor (iter = 2:d.Iterations,4)
        TanovaEffectSize(:,:,:,iter) = DoAllEffectSizesRandomized(InData,NewDesign',iDesign,SelIndFeature,DoGroup,DoF1,DoF2,bNormalize,d.ContBetween,d.ContF1,NoXing);        
         ppm.increment();
    end
else
    
    if nargin < 2
         waitbar(1/d.Iterations,h);
        set(h,'Name',sprintf('Remaining time: %01.0f:%02.0f min',floor(toc()*(d.Iterations)/60),rem(toc()*(d.Iterations-1),60)));
    else
         ShowProgress(1/d.Iterations,h);
    end

    for iter = 2:d.Iterations
        TanovaEffectSize(:,:,:,iter) = DoAllEffectSizesRandomized(InData,NewDesign',iDesign,SelIndFeature,DoGroup,DoF1,DoF2,bNormalize,d.ContBetween,d.ContF1,NoXing);        
        if nargin < 2
             waitbar(iter/d.Iterations,h);
             set(h,'Name',sprintf('Remaining time: %01.0f:%02.0f min',floor(toc()*(d.Iterations/iter-1)/60),rem(toc()*(d.Iterations/iter-1),60)));
        else
            ShowProgress(iter/d.Iterations,h);
        end
    end
end


if MeanInterval
    res.TanovaEffectSize = TanovaEffectSize;
else
    res.TanovaEffectSize = zeros(2,4,size(d.V,4),d.Iterations);
    res.TanovaEffectSize(:,:,d.StartFrame:d.EndFrame,:) = TanovaEffectSize;
end

Rank = 1:d.Iterations;
res.PTanova = ones(size(res.TanovaEffectSize)) * d.Iterations;

if d.MeanInterval
    t1 = 1;
    t2 = 1;
else
    t1 = d.StartFrame;
    t2 = d.EndFrame;
end

for i = 1:2
    for j = 1:4
        for t = t1:t2%size(d.V,4);
            if max(squeeze(res.TanovaEffectSize(i,j,t,:))) == 0
%24.2.          d.PTanova(i,j,t,:) = d.Iterations;
                res.PTanova(i,j,t,:) = d.Iterations;
            else
%24.2.          [mx,order] = sort(squeeze(d.TanovaEffectSize(i,j,t,:)),'descend');
%24.2.          d.PTanova(i,j,t,order) = Rank;
                [mx,order] = sort(squeeze(res.TanovaEffectSize(i,j,t,:)),'descend');
                res.PTanova(i,j,t,order) = Rank;
            end
        end
    end
end

%24.2   d.PTanova = d.PTanova ./ d.Iterations;
res.PTanova = res.PTanova ./ d.Iterations;
if (nargin < 2)
    close(h);
end

if d.DoGFP
    d.GFPTanovaEffectSize = res.TanovaEffectSize;
    d.GFPPTanova          = res.PTanova;
else
    d.TanovaEffectSize = res.TanovaEffectSize;
    d.PTanova          = res.PTanova;
end

result = d;

function res = DoAllEffectSizesRandomized(data,Design,iDesign,Group,DoGroup,DoF1,DoF2,Normalize,ContGroup,ContF1,NoXing)

r_in = zeros(size(data)); 
for j = 1:size(data,1)
    r_in(j,:,:,:) = data(j,PermDesign(Design,NoXing),:,:);
end
rIndFeat = Group(randperm(numel(Group)));

res = DoAllEffectSizes(r_in,Design,iDesign,rIndFeat,DoGroup,DoF1,DoF2,Normalize,ContGroup,ContF1);


function efs = DoAllEffectSizes(data,Design,iDesign,Group,DoGroup,DoF1,DoF2,Normalize,ContGroup,ContF1)

if Normalize
    data = NormDimL2(data,3);
end

if (numel(size(data)) == 3)
    DoAverage = true;
else
    DoAverage = false;
end

if DoAverage
    t2 = 1;
    efs = nan(2,4);
else
    t2 = size(data,4);
    efs = nan(2,4,size(data,4));
end

for t = 1:t2 %   size(data,4);
    if DoAverage
        in = data;
    else
        in = data(:,:,:,t);
    end
    
    %    'Take out the mean'
    [es,mmBL] = OLES(in,ones(1,size(data,2)),0);
    noBL_M = in - mmBL;
    efs(1,1,t) = es;
    
%    'Main group effect'
    if DoGroup
        if ContGroup
            [es,mmBLG] = OLESCG(noBL_M,ones(1,size(data,2)),Group,Normalize);
        else
            [es,mmBLG] = OLESG(noBL_M,ones(1,size(data,2)),Group,Normalize);
        end
        noBL_MG = noBL_M - mmBLG;
        efs(2,1,t) = es;
    end
%    'Factor 1'
    if DoF1
        if ContF1 == 0
            [es,mmBLF1] = OLES(noBL_M,Design(:,1),Normalize);
        else
            [es,mmBLF1] = OLESC(noBL_M,Design(:,1),Normalize);
        end
%        noBL_MF1 = noBL_M - mmBLF1;
        efs(1,2,t) = es;
    end
%    'Factor 1 and Group'
    if (DoF1 && DoGroup)
        if ContGroup
            if ContF1 == 0
                [es,mmBLGF1] = OLESCG(noBL_MG-mmBLF1,Design(:,1),Group,Normalize);
            else
                [es,mmBLGF1] = OLESCCG(noBL_MG-mmBLF1,Design(:,1),Group,Normalize);
            end
        else
            if ContF1 == 0
                [es,mmBLGF1] =  OLESG(noBL_MG-mmBLF1,Design(:,1),Group,Normalize);
            else
                [es,mmBLGF1] =  OLESGDC(noBL_MG-mmBLF1,Design(:,1),Group,Normalize);
            end
        end
  %      noBL_MGF1 = noBL_MG-mmBLF1-mmBLGF1;
        efs(2,2,t) = es;
    end
    % Factor 2
    if DoF2
        [es,mmBLF2] = OLES(noBL_M,Design(:,2),Normalize);
 %       noBL_MF2 = noBL_M - mmBLF2;
        efs(1,3,t) = es;
    end
    % Factor 2 and Group
    if (DoF2 && DoGroup)
        if ContGroup
            [es,mmBLGF2] = OLESCG(noBL_MG-mmBLF2,Design(:,2),Group,Normalize);
        else
            [es,mmBLGF2] = OLESG(noBL_MG-mmBLF2,Design(:,2),Group,Normalize);
        end
%        noBL_MGF2 = noBL_MG-mmBLF2-mmBLGF2;
        efs(2,3,t) = es;
    end
    % Factor 1 * Factor 2
    if (DoF1 && DoF2)
        [es,mmBLGF1F2] = OLES(noBL_M-mmBLF1-mmBLF2,iDesign,Normalize);
        efs(1,4,t) = es;
    end
    % Factor 1 * Factor 2 * Group
    if (DoF1 && DoF2 && DoGroup)
        if ContGroup
            efs(2,4,t) = OLESCG(noBL_MG-mmBLF1-mmBLF2-mmBLGF1-mmBLGF2-mmBLGF1F2,iDesign,Group,Normalize);
        else
            efs(2,4,t) =  OLESG(noBL_MG-mmBLF1-mmBLF2-mmBLGF1-mmBLGF2-mmBLGF1F2,iDesign,Group,Normalize);
        end
    end
end


function [es,mmap] = OLESCG(in,Design,Group,nFlag)

Group = normr((Group(:)-mean(Group))');

dLevel = unique(Design);
ndLevel = numel(dLevel);

cm = zeros(ndLevel,size(in,3));

for l = 1:ndLevel
    in_s = in(:,Design == dLevel(l),:);
    m = squeeze(mean(in_s,2));
    cm(l,:) = Group *m;
end

if nFlag
    disp('Dissimilarity not inplemented for continuous predictors');
end

% Normalization not implemented

es = sqrt(cm.*cm);
es = sum(es(:));

cmp = zeros(size(in,1),ndLevel,size(in,3));
parfor s = 1:size(in,1)
    cmp(s,:,:) = cm * Group(s);
end

if nargout > 1
    mmap = zeros(size(in));
    for l = 1:ndLevel
        idx = find(Design == dLevel(l));
        mmap(:,idx,:) = repmat(cmp(:,l,:),[1 numel(idx),1]);
    end
end

function [es,mmap] = OLESG(in,Design,Group,nFlag)
%nSub = size(in,1);
Level = unique(Group);
nLevel = numel(Level);
dLevel = unique(Design);
ndLevel = numel(dLevel);

LMap = zeros(nLevel,ndLevel,size(in,3));
for l = 1:nLevel
    for  l2 = 1:ndLevel
        in_s = in(Group == Level(l),Design == dLevel(l2),:);
        nCObs = size(in_s,1) * size(in_s,2);
        maps = reshape(in_s,nCObs,size(in_s,3));
        LMap(l,l2,:) = mean(maps); 
    end
end

if nFlag == 1
    LMap2 = NormDimL2(LMap,3);
    es = sqrt(sum(LMap2.*LMap2,3));
    es = sum(es(:));
else
    es = sqrt(sum(LMap.*LMap,3));
    es = sum(es(:));
end

if nargout > 1
    mmap = zeros(size(in));
    for l = 1:nLevel
        idx = find(Group == Level(l));
        for l2 = 1:ndLevel
            idx2 = find(Design == dLevel(l2));
            mmap(idx,idx2,:) = repmat(LMap(l,l2,:),[numel(idx),numel(idx2),1]);
        end
    end
end



function [es,mmap] = OLES(in,Design,nFlag)
nSub = size(in,1);
Level = unique(Design);
nLevel = numel(Level);

LMap = zeros(1,nLevel,size(in,3));

for l = 1:nLevel
    in_s = in(:,Design == Level(l),:);
    nCObs = size(in_s,1) * size(in_s,2);
    maps = reshape(in_s,nCObs,size(in_s,3));
    LMap(1,l,:) = mean(maps); 
end

if nFlag == 1
    LMap2 = NormDimL2(LMap,3);
    es = sum((sqrt(sum(LMap2.*LMap2,3))));
else
    es = sum((sqrt(sum(LMap.*LMap,3))));
end

if nargout > 1
    mmap = zeros(size(in));
    for l = 1:nLevel
        idx = find(Design == Level(l));
        mmap(:,idx,:) = repmat(LMap(1,l,:),[nSub,numel(idx),1]);        
    end
end

function [es,mmap] = OLESC(in,Design,nFlag)

in_s = squeeze(mean(in,1));
Design = normr((Design(:)-mean(Design))');

LMap(1,:,:) = Design * in_s;
if nFlag
    disp('Dissimilarity not inplemented for continuous predictors');
end

es = sqrt(sum(LMap.*LMap,3));

mmap = repmat(LMap,[size(in,1) size(in,2) 1]);


function [es,mmap] = OLESCCG(in,Design,Group,nFlag)

Group = normr((Group(:)-mean(Group)));
Design = normr((Design(:)-mean(Design))');

GroupR  = repmat(Group ,[1 size(in,2) size(in,3)]);
DesignR = repmat(Design,[size(in,1),1,size(in,3)]); 

TotDesignR = DesignR .* GroupR;

LMap = mean(mean(in.*TotDesignR,1),2);

if nFlag
    disp('Dissimilarity not inplemented for continuous predictors');
end

% Normalization not implemented

es = sqrt(LMap.*LMap);
es = sum(es(:));

if nargout > 1
    LMapR = repmat(LMap,[size(in,1) size(in,2) 1]);
    mmap = LMapR .*TotDesignR;
end


function [es,mmap] = OLESGDC(in,Design,Group,nFlag)
%nSub = size(in,1);
Level = unique(Group);
nLevel = numel(Level);

Design = normr((Design(:)-mean(Design))');

LMap = zeros(nLevel,size(in,3));
for l = 1:nLevel
    in_s = squeeze(mean(in(Group == Level(l),:,:)));
    LMap(l,:) = Design * in_s;
end

if nFlag
    disp('Dissimilarity not inplemented for continuous predictors');
end

es = sqrt(LMap.*LMap);
es = sum(es(:));

if nargout > 1
    mmap = zeros(size(in));
    for l = 1:nLevel
        idx = find(Group == Level(l));
        dmap(1,:,:) = Design' * squeeze(LMap(l,:));
        mmap(idx,:,:) = repmat(dmap,[numel(idx),1,1]);
    end
end

