% K-means Clustering
% (as described in Cartool)
%
% Result = lab_kmeans(data,settings)
%
%   data                 = nchans x ntimeframes
%   settings.maxclusters = maximal number of clusters
%                          (default = 20)
%   settings.MaxIter     = number of iterations per number of clusters
%
% written by F. Hatz 2014

function Result = lab_kmeans(data,settings)

disp('   K-means clustering')

if ~isfield(settings,'minclusters')
    settings.minclusters = 2;
end
if ~isfield(settings,'maxclusters')
    settings.maxclusters = 20;
end
if ~isfield(settings,'MaxIter')
    settings.MaxIter = 300;
end
if ~isfield(settings,'maxerror')
    settings.maxerror = 1e-9;
end

Ntr = size(data,2);
Nchans = size(data,1);
NumDots = floor(settings.MaxIter/20);
GFP2 = sum(data.^2,1) / Nchans;

ClusterNr = settings.minclusters:settings.maxclusters;
numclusters = length(ClusterNr);
Result.Nr = ClusterNr;
Result.Clusters = zeros(numclusters,Ntr);
Result.GEV = zeros(1,numclusters);
Result.SumError = ones(1,numclusters);
Result.Niter = zeros(1,numclusters);
Result.Distance = zeros(1,numclusters);
Result.CORR = zeros(numclusters,Ntr);
for Nclust = 1:numclusters
    fprintf(['     calculate K-means (' num2str(ClusterNr(Nclust)) ' Clusters)'])
    % initialize
    TemplateOLD = zeros(Nchans,ClusterNr(Nclust));
    Template = data(:,randperm(Ntr,ClusterNr(Nclust)));
    CORR = corr(Template,data);
    [~,IDX] = max(CORR.^2,[],1);
    for j = 1:ClusterNr(Nclust)
        Template(:,j) = mean(data(:,IDX==j) .* repmat(sign(CORR(j,IDX==j)),Nchans,1),2);
    end
    % iterate for best result
    for Niter = 1:settings.MaxIter
        % fit data to template
        CORR = corr(Template,data);
        [EV,IDX] = max(CORR.^2,[],1);
        GEV = sum(EV .* GFP2) / sum(GFP2);
        if GEV > Result.GEV(1,Nclust)
            Result.Clusters(Nclust,:) = IDX;
            Result.CORR(Nclust,:) = EV;
            Result.CorrAll{Nclust} = CORR.^2;
            Result.Template{Nclust} = Template;
            Result.GEV(1,Nclust) = GEV;
            Result.Distance(1,Nclust) = sum((abs(Template(:)-TemplateOLD(:))).^2).^0.5;
            Result.Niter(1,Nclust) = Niter;
        end
        % check max error
        if Result.Distance(1,Nclust) < settings.maxerror
            break
        end
        % calculate new template
        TemplateOLD = Template;
        for j = 1:ClusterNr(Nclust)
            Template(:,j) = mean(data(:,IDX==j) .* repmat(sign(CORR(j,IDX==j)),Nchans,1),2);
        end
        if mod(Niter,NumDots) == 0
            fprintf('.');
        end
    end
    disp(':');
end

end
