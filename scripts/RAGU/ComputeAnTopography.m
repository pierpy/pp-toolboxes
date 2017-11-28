function res = ComputeAnTopography(out,h)

% Copyright 2009-2011 Thomas Koenig
% distributed under the terms of the GNU AFFERO General Public License

if nargin < 2
    h = waitbar(0,'Computing, please wait...');
end

if ~out.ContBetween
    b = unique(out.IndFeature(~isnan(out.IndFeature)));
    for i = 1:numel(b)
        Group{i} = find(out.IndFeature == b(i));
    end
else
    Group{1} = 1:numel(out.IndFeature(~isnan(out.IndFeature)));
end 


if out.MeanInterval == 1
   DataToUse = mean(out.V(:,:,:,out.StartFrame:out.EndFrame),4);
   StartPoint = 1;
   EndPoint   = 1;
   nPoints    = 1;
else
    DataToUse = out.V(:,:,:,out.StartFrame:out.EndFrame);
    StartPoint = 1;
    EndPoint   = out.EndFrame - out.StartFrame+1;
    nPoints    = EndPoint - StartPoint + 1;
end

if out.DoAnTopFFT == 0
    out.pTopCons    = ones(numel(Group),size(out.V,2),nPoints);
    pAllTopCons = ones(numel(Group),size(out.V,2),nPoints,out.Iterations);
    out.MeanGFP     = ones(numel(Group),size(out.V,2),nPoints);
    
    tic
    totTests = (nPoints*size(out.V,2));
    for c = 1:size(out.V,2)
        for t = StartPoint:EndPoint
            for grp = 1:numel(Group)
                in = squeeze(out.V(Group{grp},c,:,(t+out.StartFrame-1)));
                [gt,pt,AllP] = AnTop(in,out.Iterations,'q');
                pAllTopCons(grp,c,t,:) = AllP;
                out.pTopCons(grp,c,t) = pt;
                out.MeanGFP(grp,c,t)  = gt;

            end
            if nargin < 2
                n = ((c-1) * nPoints + t);
                waitbar(n/totTests,h);
                set(h,'Name',sprintf('Remaining time: %01.0f:%02.0f min',floor(toc()*(totTests/n-1)/60),rem(toc()*(totTests/n-1),60)));
            else
                ShowProgress(((c-1) * nPoints + t)/(nPoints*size(out.V,2)),h);
            end
        end
    end
else
    VFFT = fft(out.V,[],4);
    nPts = floor(size(VFFT,4) / 2+1);
    VFFT = VFFT(:,:,:,1:nPts);
    out.pTopCons    = ones(numel(Group),size(VFFT,2),size(VFFT,4));
    pAllTopCons = ones(numel(Group),size(out.V,2),size(out.V,4),out.Iterations);
    out.MeanGFP     = ones(numel(Group),size(VFFT,2),size(VFFT,4));
    tic
    totTests = (size(VFFT,4)*size(VFFT,2));
    for c = 1:size(VFFT,2)
        for t = 1:size(VFFT,4)
            for grp = 1:numel(Group)
                if nargin < 2
                    n = ((c-1) * size(VFFT,4) + t);
                    waitbar(n/totTests,h);
                    set(h,'Name',sprintf('Remaining time: %01.0f:%02.0f min',floor(toc()*(totTests/n-1)/60),rem(toc()*(totTests/n-1),60)));
                else
                    ShowProgress(((c-1) * size(VFFT,4) + t)/(size(VFFT,4)*size(VFFT,2)),h);
                end
                in = squeeze(VFFT(Group{grp},c,:,t));
                [gt,pt,AllP] = AnTop(in,out.Iterations,'q');
                pAllTopCons(grp,c,t,:) = AllP;
                out.pTopCons(grp,c,t) = pt;
                out.MeanGFP(grp,c,t)  = gt;
            end
        end
    end
end 
 

for g = 1:size(out.pTopCons,1)
    for c = 1:size(out.pTopCons,2)
        hit = squeeze(pAllTopCons(g,c,:,:)) < out.Threshold;
        hitCount = squeeze(sum(hit,1));
        out.TCTPHitCount(g,c) = sum((hitCount >= hitCount(1))) / size(hit,2);
        out.TCTHitDuration{g,c} = RaguClustSize(squeeze(hit(:,:)));
        p = 1-(cumsum(out.TCTHitDuration{g,c})/sum(out.TCTHitDuration{g,c}));
        crit = find(p < out.Threshold);
        if ~isempty(crit)
            out.TCTCritDuration(g,c) = crit(1);
        else
            out.TCTCritDuration(g,c) = size(out.V,4);
        end
    end
end

if nargin < 2
    close(h);
end

res = out;