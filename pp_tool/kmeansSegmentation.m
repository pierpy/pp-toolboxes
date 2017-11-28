clear 
close all
clc

files_path='/Users/pp/Documents/MATLAB/microstates_pp/dati/Rest9.DECI.DECIMATE_MAX.ep';



data = load(files_path);
data = data';
Ntr = size(data,2);
Nchans = size(data,1);
 
%average reference

%

settings.maxclusters = 10; 
settings.minclusters = 2;
settings.MaxIter = 300;
settings.maxerror = 1e-9;


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
    CORR2 = CORR.^2;
    [~,IDX] = max(CORR2,[],1);
    for j = 1:ClusterNr(Nclust)
        Template(:,j) = mean(data(:,IDX==j) .* repmat(sign(CORR(j,IDX==j)),Nchans,1),2);
    end
    % iterate for best result
    for Niter = 1:settings.MaxIter
        % fit data to template
        CORR = corr(Template,data);
        CORR2 = CORR.^2;
        [EV,IDX] = max(CORR2,[],1);
        GEV = sum(EV .* GFP2) / sum(GFP2);
        cv_coeff = (Nchans-1)/(Nchans-1-ClusterNr(Nclust));
        cv_coeff = cv_coeff^2;
        if GEV > Result.GEV(1,Nclust)
            Result.Clusters(Nclust,:) = IDX;
            Result.CORR(Nclust,:) = EV;
            Result.CorrAll{Nclust} = CORR.^2;
            Result.Template{Nclust} = Template;
            Result.GEV(1,Nclust) = GEV;
            Result.Distance(1,Nclust) = sum((abs(Template(:)-TemplateOLD(:))).^2).^0.5;
            Result.Niter(1,Nclust) = Niter;
            %% calcolo CV
            templateTemp = Result.Template{1,Nclust};
            templateTemp2 = templateTemp;
            %normalizzo i template
            for j = 1:size(templateTemp,2) 
                templateTemp2(:,j) = templateTemp(:,j) ./ (sum(templateTemp(:,j).^2)).^0.5;
            end
            datatmp = data(:,Result.Clusters(Nclust,:)>0);
            Clusterstmp = Result.Clusters(Nclust,Result.Clusters(Nclust,:)>0);
            tmp = zeros(1,size(datatmp,2));
            for j = 1:size(datatmp,2)
                 tmp(1,j) = datatmp(:,j)'*datatmp(:,j) - (templateTemp2(:,Clusterstmp(1,j))'*datatmp(:,j))^2;
            end
            Result.CV(1,Nclust) = (sum(tmp)/(Ntr*(Nchans-1)))*cv_coeff;
            %% calcolo W (dispersion measure)
            for j = 1:ClusterNr(Nclust)

                D(j)= sum(sum(triu(squareform(pdist(data(:,Clusterstmp==j),'euclidean')))))/(2*size(data(:,Clusterstmp==j),2));
                 %D(j)= sum(sum(squareform(pdist(datatmp(:,Clusterstmp==j)))))/size(datatmp(:,Clusterstmp==j),2);
            end
%             Wq = (1/2)*sum(D);
%             Result.Wq(1,Nclust)=Wq;
                Result.D{1,Nclust}=D;

                Result.W(1,Nclust)=sum(D);
                Result.Mq(1,Nclust) = sum(D)*(ClusterNr(Nclust)^(2/Nchans));
%                 Result.W2(1,Nclust)=sum(D2);
                clearvars D
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

M=Result.Mq;
for i =1:size(M,2)-1
    d(i)=M(i)-M(i+1);
end
for jj=2:size(d,2)
    if d(jj-1)<0 | d(jj-1)<d(jj)
        KL(jj)=0;
    else
        KL(jj)=(d(jj-1)-d(jj))/M(jj-1);
    end
    
end