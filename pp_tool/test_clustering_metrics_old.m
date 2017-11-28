clear 
close all
clc

%[eeg cp idx] = generateData(1, 0.5, 4, 15, 15, 5, 1, 2, 3000);

% DataSet = [5*randn(1,3000);20+3*randn(1,3000);120+25*randn(1,3000)];
files_path='/Users/pp/Documents/MATLAB/microstates_pp/dati/gfpMax_Restcermat_microstati.ep';

% 
% 
% 
data = load(files_path);
pmode=1;
minclusters=2;
maxclusters=10;
ClusterNr = minclusters:maxclusters;
numclusters = length(ClusterNr);
[Ntr, Nchans]=size(data);
for q=1:numclusters
[maps,ind,assigned_cluster_ampl,exp_var] = NMicrostates(data,ClusterNr(q),pmode,20);
            %% calcolo CV
%             cv_coeff = (Nchans-1)/(Nchans-1-ClusterNr(q));
%             cv_coeff = cv_coeff^2;
%             templateTemp = maps;
             Clusterstmp = ind;
%             tmp = zeros(size(data,1),1);
%             for j = 1:size(data,1)
%                  tmp(j,1) = data(j,:)*data(j,:)' - (templateTemp(Clusterstmp(j,1),:)*data(j,:)').^2;
%             end
%             Result.CV(q) = (sum(tmp)/(Ntr*(Nchans-1)))*cv_coeff;
            [Result.CV(q)] = CV(data,maps,ind,ClusterNr(q));
            
            %% calcolo W (dispersion measure)
            for j = 1:ClusterNr(q)

                D(j)= sum(sum(triu(squareform(pdist(data(Clusterstmp==j,:),'euclidean')))))/(2*size(data(Clusterstmp==j,:),2));
                
            end
            Result.D{q}=D;
            Result.W(q)=sum(D); 
             [Result.W1(q)] = W(data,ind,ClusterNr(q));
            clearvars D
            % result.template{n_mod}=maps;
            % result.clusters{n_mod}=ind;
            % result.amplitude{n_mod}=assigned_cluster_ampl;
            % result.exp_var=exp_var;
            clearvars maps ind assigned_cluster_ampl exp_var
end
%data = data';
% EEG.chanlocs = readlocs('coord.xyz');
% figure;
% 
% for k=1:4
% subplot(1,4,k);
% topoplot(result.template{1, 4}(k,:), EEG.chanlocs)
% title(num2str(k))
% end
% data=DataSet';
%% N-microstates algorithm
% data=data';
% 
% settings.maxclusters = 10; 
% settings.minclusters = 2;
% settings.MaxIter = 300;
% settings.maxerror = 1e-9;
% 
% 
% Ntr = size(data,2);
% Nchans = size(data,1);
% NumDots = floor(settings.MaxIter/20);
% GFP2 = sum(data.^2,1) / Nchans;
% 
% ClusterNr = settings.minclusters:settings.maxclusters;
% numclusters = length(ClusterNr);
% Result.Nr = ClusterNr;
% Result.Clusters = zeros(numclusters,Ntr);
% Result.GEV = zeros(1,numclusters);
% Result.SumError = ones(1,numclusters);
% Result.Niter = zeros(1,numclusters);
% Result.Distance = zeros(1,numclusters);
% Result.CORR = zeros(numclusters,Ntr);
% 
% for Nclust = 1:numclusters
%     fprintf(['     calculate K-means (' num2str(ClusterNr(Nclust)) ' Clusters)'])
%     % initialize
%     TemplateOLD = zeros(Nchans,ClusterNr(Nclust));
%     Template = data(:,randperm(Ntr,ClusterNr(Nclust)));
%     CORR = corr(Template,data);
%     CORR2 = CORR.^2;
%     [~,IDX] = max(CORR2,[],1);
%     for j = 1:ClusterNr(Nclust)
%         Template(:,j) = mean(data(:,IDX==j) .* repmat(sign(CORR(j,IDX==j)),Nchans,1),2);
%     end
%     % iterate for best result
%     for Niter = 1:settings.MaxIter
%         % fit data to template
%         CORR = corr(Template,data);
%         CORR2 = CORR.^2;
%         [EV,IDX] = max(CORR2,[],1);
%         GEV = sum(EV .* GFP2) / sum(GFP2);
%         cv_coeff = (Nchans-1)/(Nchans-1-ClusterNr(Nclust));
%         cv_coeff = cv_coeff^2;
%         if GEV > Result.GEV(1,Nclust)
%             Result.Clusters(Nclust,:) = IDX;
%             Result.CORR(Nclust,:) = EV;
%             Result.CorrAll{Nclust} = CORR.^2;
%             Result.Template{Nclust} = Template;
%             Result.GEV(1,Nclust) = GEV;
%             Result.Distance(1,Nclust) = sum((abs(Template(:)-TemplateOLD(:))).^2).^0.5;
%             Result.Niter(1,Nclust) = Niter;
%             %  calcolo CV
%             templateTemp = Result.Template{1,Nclust};
%             templateTemp2 = templateTemp;
%             %normalizzo i template
%             for j = 1:size(templateTemp,2) 
%                 templateTemp2(:,j) = templateTemp(:,j) ./ (sum(templateTemp(:,j).^2)).^0.5;
%             end
%             datatmp = data(:,Result.Clusters(Nclust,:)>0);
%             Clusterstmp = Result.Clusters(Nclust,Result.Clusters(Nclust,:)>0);
%             tmp = zeros(1,size(datatmp,2));
%             for j = 1:size(datatmp,2)
%                  tmp(1,j) = datatmp(:,j)'*datatmp(:,j) - (templateTemp2(:,Clusterstmp(1,j))'*datatmp(:,j))^2;
%             end
%             Result.CV(1,Nclust) = (sum(tmp)/(Ntr*(Nchans-1)))*cv_coeff;
%             % calcolo W (dispersion measure)
%             for j = 1:ClusterNr(Nclust)
% %                 maps_clu_j=datatmp(:,Clusterstmp==j);
% %                 nr=size(maps_clu_j,2);
% %                 dd=pdist(maps_clu_j);
% %                 dd=squareform(dd);
% %                 dd=triu(dd);
% %                 dd=sum(sum(dd));
% %                 D(j)=(1/nr)*dd;
%                 D(j)= sum(sum(triu(squareform(pdist(datatmp(:,Clusterstmp==j))))))/(2*size(datatmp(:,Clusterstmp==j),2));
% %                 D2(j)= sum(sum(squareform(pdist(datatmp(:,Clusterstmp==j)))))/size(datatmp(:,Clusterstmp==j),2);
%             end
% %             Wq = (1/2)*sum(D);
% %             Result.Wq(1,Nclust)=Wq;
% 
%                 Result.W(1,Nclust)=sum(D);
% %                 Result.W2(1,Nclust)=sum(D2);
%                 clearvars D
%         end
%         % check max error
%         if Result.Distance(1,Nclust) < settings.maxerror
%             break
%         end
%         % calculate new template
%         TemplateOLD = Template;
%         for j = 1:ClusterNr(Nclust)
%             Template(:,j) = mean(data(:,IDX==j) .* repmat(sign(CORR(j,IDX==j)),Nchans,1),2);
%         end
%         if mod(Niter,NumDots) == 0
%             fprintf('.');
%         end
%     end
%     disp(':');
% end
% figure
% plot(Result.CV)
% figure
% plot(Result.W)
%Result = lab_kmeans(segment,settings);
% Result = mskmeans(data,settings);

% [Result,Stats] = lab_clustering_stats(Result,segment);

%% simple k-means algorithm

% Nchans=size(data,2);
% Ntr=size(data,1);
% GFP = sum(data.^2,2) / Nchans;
% maxclusters = 10;
% 
% for q = 1:maxclusters
%     [clu, template]=kmeans(data,q);
%     templates{q}=template;
%     clusters(q,:)=clu;
%     %% calcolo CV
%             cv_coeff = (Nchans-1)/(Nchans-1-q);
%             cv_coeff = cv_coeff^2;
%             templateTemp = templates{q};
%             templateTemp2 = templateTemp;
%             %normalizzo i template
%             for j = 1:size(templateTemp,1) 
%                 templateTemp2(j,:) = templateTemp(j,:) ./ (sum(templateTemp(j,:).^2)).^0.5;
%             end
%             Clusterstmp = clusters(q,:);
%             tmp = zeros(size(data,1),1);
%             for j = 1:size(data,1)
%                  tmp(j,1) = data(j,:)*data(j,:)' - (templateTemp2(Clusterstmp(1,j),:)*data(j,:)').^2;
%             end
%             Result.CV(q) = (sum(tmp)/(Ntr*(Nchans-1)))*cv_coeff;
%             
%            % calcolo W (dispersion measure)
%             for j = 1:q
% 
%                 D(j)= sum(sum(triu(squareform(pdist(data(Clusterstmp==j,:),'euclidean')))))/(2*size(data(Clusterstmp==j,:),2));
%                 
%             end
%             
%                 Result.D{q}=D;
%                 Result.W(q)=sum(D);
%                 Result.Mq(q) = sum(D)*(q^(2/Nchans));
% 
%                 clearvars D
% 
% end
% Mq=Result.Mq;
% for ii = 1: size(Mq,2)-1
%     d(ii) = Mq(ii) - Mq(ii+1);
% end
% 
% for jj=2:size(d,2)
% %     if d(jj-1)<0 | d(jj-1)<d(jj)
% %         KL(jj)=0;
% %     else
%         KL(jj)=(d(jj-1)-d(jj))/Mq(jj-1);
% %     end
%     
% end
% plot(Result.W)
% figure
% plot(KL)