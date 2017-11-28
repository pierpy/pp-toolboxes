function [Result,Stats] = lab_clustering_stats(Result,data)
    
if ~isfield(Result,'Template') | ~isfield(Result,'Clusters') | size(Result.Clusters,2) ~= size(data,2)
    Stats = [];
    return
end

ClusterNr = Result.Nr;
numclusters = length(ClusterNr);
nChans = size(data,1);
GFP2 = sum(data.^2,1) / nChans;
Stats = [];
fprintf('   Calculate Microstates-Stats');
for Nclust = 1:numclusters
    Stats(Nclust).Nr = ClusterNr(Nclust); %#ok<AGROW>
    Stats(Nclust).CORR = zeros(1,size(data,2)); %#ok<AGROW>
    Stats(Nclust).ClusterTopo = zeros(nChans,ClusterNr(Nclust)); %#ok<AGROW>
    Clusters = Result.Clusters(Nclust,:);
    for j = 1:ClusterNr(Nclust)
        Idx = find(Clusters==j);
        if ~isempty(Idx)
            datatmp = data(:,Idx);
            CORR = corr(Result.Template{Nclust}(:,j),datatmp);
            Stats(Nclust).CORR(1,Idx) = CORR; %#ok<AGROW>
            Stats(Nclust).MeanCorr(1,j) = mean(abs(CORR)); %#ok<AGROW>
            Stats(Nclust).MeanEV(1,j) = mean(CORR.^2); %#ok<AGROW>
            [Stats(Nclust).MaxCorr(1,j),tmp] = max(abs(CORR)); %#ok<AGROW>
            tmp = Idx(tmp);
            Stats(Nclust).MaxCorrIdx(1,j) = tmp(1); %#ok<AGROW>
            Stats(Nclust).ClusterGEV(1,j) = sum(CORR.^2 .* GFP2(Idx)) / sum(GFP2(Idx)); %#ok<AGROW>
            
            tcorr = sign(corr(mean(datatmp,2),datatmp));
            datatmp = datatmp .* repmat(tcorr(:)',size(datatmp,1),1);
            Stats(Nclust).ClusterTopo(:,j) = mean(datatmp,2); %#ok<AGROW>
            [Stats(Nclust).MaxGfp(1,j),tmp] = max(GFP2(Idx)); %#ok<AGROW>
            tmp = Idx(tmp);
            Stats(Nclust).MaxGfpIdx(1,j) = tmp(1); %#ok<AGROW>
        else
            Stats(Nclust).MeanCorr(1,j) = 0; %#ok<AGROW>
            Stats(Nclust).MeanEV(1,j) = 0; %#ok<AGROW>
            Stats(Nclust).MaxCorr(1,j) = 0; %#ok<AGROW>
            Stats(Nclust).MaxCorrIdx(1,j) = 0; %#ok<AGROW>
            Stats(Nclust).ClusterGEV(1,j) = 0; %#ok<AGROW>
            Stats(Nclust).ClusterTopo(:,j) = 0; %#ok<AGROW>
            Stats(Nclust).MaxGfp(1,j) = 0; %#ok<AGROW>
            Stats(Nclust).MaxGfpIdx(1,j) = 0; %#ok<AGROW>
        end
    end
    Stats(Nclust).GEV = sum(Stats(Nclust).CORR.^2 .* GFP2) / sum(GFP2); %#ok<AGROW>
    Stats(Nclust).KL = (ClusterNr(Nclust)^(2/size(data,2))) * mean(Stats(Nclust).CORR.^2); %#ok<AGROW>
    Stats(Nclust).MeanVariance = mean(Stats(Nclust).CORR.^2); %#ok<AGROW>
    Result.MeanVariance(1,Nclust) = Stats(Nclust).MeanVariance;
    Result.GEV(1,Nclust) = Stats(Nclust).GEV;
    Result.KL(1,Nclust) = Stats(Nclust).KL;
    
    % calculate cross validation
    Template = Result.Template{Nclust};
    Template2 = Template;
    for j = 1:size(Template,2)
        Template2(:,j) = Template(:,j) ./ (sum(Template(:,j).^2)).^0.5;
    end
    datatmp = data(:,Result.Clusters(Nclust,:)>0);
    Clusterstmp = Result.Clusters(Nclust,Result.Clusters(Nclust,:)>0);
    tmp = zeros(1,size(datatmp,2));
    for j = 1:size(datatmp,2)
        tmp(1,j) = datatmp(:,j)'*datatmp(:,j) - (Template2(:,Clusterstmp(1,j))'*datatmp(:,j))^2;
    end
    Result.CV(1,Nclust) = sum(tmp) / (size(datatmp,2)*(size(datatmp,1)-1));
    Result.CV(1,Nclust) = Result.CV(1,Nclust) * ((size(datatmp,1)-1)/(size(datatmp,1)-1-ClusterNr(Nclust)))^2;
    Stats(Nclust).CV = Result.CV(1,Nclust); %#ok<AGROW>
    clearvars datatmp Template2 Clusterstmp
    
    % calculate cluster stats
    tmp = find(abs(diff(Clusters))>0);
    tmp = tmp(:)';
    EndMicro = [tmp size(data,2)];
    StartMicro = [1 tmp+1];
    clearvars tmp
    IndexMicro = Clusters(1,StartMicro);
    tmp = [0 Clusters];
    PreviousMicro = tmp(1,StartMicro);
    tmp = [Clusters 0];
    NextMicro = tmp(1,EndMicro+1);
    clearvars tmp
    for j = 1:ClusterNr(Nclust)
        StartMicro2 = StartMicro(IndexMicro==j);
        EndMicro2 = EndMicro(IndexMicro==j);
        if ~isempty(StartMicro2)
            TF = EndMicro2 - StartMicro2 +1;
            Stats(Nclust).MeanDuration(1,j) = mean(TF); %#ok<AGROW>
            Stats(Nclust).MinDuration(1,j) = min(TF); %#ok<AGROW>
            Stats(Nclust).MaxDuration(1,j) = max(TF); %#ok<AGROW>
            Stats(Nclust).StdDuration(1,j) = std(TF); %#ok<AGROW>
            Stats(Nclust).SumDuration(1,j) = sum(TF); %#ok<AGROW>
            PreviousMicro2 = PreviousMicro(IndexMicro==j);
            NextMicro2 = NextMicro(IndexMicro==j);
            for i = 1:ClusterNr(Nclust)
                Stats(Nclust).PreviousCluster(i,j) = sum(PreviousMicro2 == i); %#ok<AGROW>
                Stats(Nclust).NextCluster(i,j) = sum(NextMicro2 == i); %#ok<AGROW>
            end
        else
            Stats(Nclust).MeanDuration(1,j) = 0; %#ok<AGROW>
            Stats(Nclust).MinDuration(1,j) = 0; %#ok<AGROW>
            Stats(Nclust).MaxDuration(1,j) = 0; %#ok<AGROW>
            Stats(Nclust).StdDuration(1,j) = 0; %#ok<AGROW>
            Stats(Nclust).SumDuration(1,j) = 0; %#ok<AGROW>
            for i = 1:ClusterNr(Nclust)
                Stats(Nclust).PreviousCluster(i,j) = 0; %#ok<AGROW>
                Stats(Nclust).NextCluster(i,j) = 0; %#ok<AGROW>
            end
        end
    end
    fprintf('.');
end
disp(':');

end