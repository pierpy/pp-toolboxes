% Topographic Atomize & Agglomerate Hierarchical Clustering (T-AAHC)
% (as described in Cartool)
%
% Result = lab_taahc(data,settings)
%
%   data                 = nchans x ntimeframes
%   settings.maxclusters = maximal number of clusters to store
%                                 (default = 20)
%
% written by F. Hatz 2014

function Result = lab_taahc(data,settings)

disp('   T-AAHC clustering')

if ~exist('settings','var') | ~isfield(settings,'minclusters')
    settings.minclusters = 2;
end
if ~exist('settings','var') | ~isfield(settings,'maxclusters')
    settings.maxclusters = 20;
end

Ntr = size(data,2);
Nchans = size(data,1);
NumDots = floor(Ntr/20);
Clusters = 1:Ntr;
Nclust = max(Clusters);
MeanEV = ones(1,Ntr);
Template = data;
Tnr = ones(1,Ntr);
GFP2 = sum(data.^2,1) / Nchans;

fprintf('     calculate T-AAHC')
Result.Nr = [];
Result.Clusters = [];
Result.Template = {};
while Nclust > 1
    if Nclust <= settings.maxclusters
        Result.Nr = [Nclust Result.Nr];
        Result.Clusters = cat(1,Clusters,Result.Clusters);
        Result.Template= [{Template};Result.Template];
    end
    if Nclust > settings.minclusters
        tmp = find(MeanEV == min(MeanEV),1,'last');
        Cindex = setdiff(1:Nclust,tmp);
        Template = Template(:,Cindex);
        Tnr = Tnr(:,Cindex);
        MeanEV = MeanEV(:,Cindex);
        tmp = find(Clusters == tmp);
        for i = tmp
            CORR = corr(data(:,i),Template);
            Tidx = find(CORR.^2 == max(CORR.^2),1,'first');
            Clusters(i) = Cindex(Tidx);
            Template(:,Tidx) = (Tnr(1,Tidx)*Template(:,Tidx) + sign(CORR(Tidx))*data(:,i)) / (Tnr(1,Tidx) + 1);
            Tnr(1,Tidx) = Tnr(1,Tidx) + 1;
            MeanEV(1,Tidx) = mean(corr(Template(:,Tidx),data(:,Clusters==Cindex(Tidx))).^2);
        end
        [~,~,Clusters] = unique(Clusters);
        Nclust = max(Clusters);
    else
        break;
    end
    if mod(Nclust,NumDots) == 0
        fprintf('.');
    end
end
disp(':');

end
