clear
close all
clc

files_path='/Users/pp/Documents/MATLAB/microstates_pp/dati/gfpMax_Restdegmat_microstati.ep';

% cd(files_path);
% files=dir;
% files(1:2)=[];


eeg = load(files_path);

%% simple k-means
nch=size(eeg,2);
nsamples=size(eeg,1);
GFP = sum(eeg.^2,2) / nch;
maxclusters = 10;

for q = 1:maxclusters
    [clu, template]=kmeans(eeg,q);
    templates{q}=template;
    clusters(q,:)=clu;
    coeff = (nch-1)/(nch-1-q);
    coeff = coeff^2;
        for kk=1:nsamples
            a(kk)=GFP(kk)*corr2(eeg(kk,:),templates{q}(clusters(q,kk),:));
            b(kk)=GFP(kk)^2;
            c(kk)=sum(eeg(kk,:).^2)-dot(templates{q}(clusters(q,kk),:),eeg(kk,:));
        end
        a = a.^2;
        gev(q)=sum(a)/sum(b);
        d(q)=sum(c)/(nsamples*(nch-1));
        cv(q)=d(q)*coeff;
        
       for i = 1:q 
           temp=pdist(eeg(find(clusters(q,:)==i),:));
           temp=squareform(temp);
           temp=triu(temp);
           temp=sum(temp);
           Dr(i)=sum(temp);
       end
       Wq(q)=(1/2*q)*sum(Dr);
       Mq(q)=Wq(q)*(q^(2/nch));
end

for ii = 1: size(Mq,2)-1
    d(ii) = Mq(ii) - Mq(ii+1);
end

for jj=2:size(d,2)
    if d(jj-1)<0 | d(jj-1)<d(jj)
        KL(jj)=0;
    else
        KL(jj)=(d(jj-1)-d(jj))/Mq(jj-1);
    end
    
end




%% cartool
eeg_c = eeg';
nch_c = size(eeg_c,1);
nsamples_c = size(eeg_c,2);
GFP_c = sum(eeg_c.^2,1) / nch_c;

settings.maxclusters = 10; 
settings.minclusters = 1;
settings.MaxIter = 300;
result = lab_kmeans(eeg_c,settings);

numCluster = result.Nr;
cluster = result.Clusters;
corr = result.CORR;
templates_c = result.Template;
%%%gev
for q = 1: size(numCluster,2)
        for kk=1:nsamples
            a_c(kk)=GFP_c(kk)*corr2(eeg_c(:,kk),templates_c{q}(:,cluster(q,kk)));
            b_c(kk)=GFP_c(kk)^2;
            c_c(kk)=sum(eeg_c(:,kk).^2)-dot(templates_c{q}(:,cluster(q,kk)),eeg_c(:,kk));
        end
        a_c = a_c.^2;
        gev_c(q)=sum(a_c)/sum(b_c);
        d_c(q)=sum(c_c)/(nsamples*(nch-1));
        cv_c(q)=d_c(q)*coeff;
        
       for i = 1:q 
           temp_c=pdist(eeg_c(:,find(clusters(q,:)==i)));
           temp_c=squareform(temp_c);
           temp_c=triu(temp_c);
           temp_c=sum(temp_c);
           Dr_c(i)=sum(temp_c);
       end
       Wq_c(q)=(1/2*q)*sum(Dr_c);
       Mq_c(q)=Wq_c(q)*(q^(2/nch));  
end
%%%
for ii = 1: size(Mq,2)-1
    d_c(ii) = Mq_c(ii) - Mq_c(ii+1);
end

for jj=2:size(d,2)
    %if d_c(jj-1)<0 | d_c(jj-1)<d_c(jj)
        %KL_c(jj)=0;
    %else
        KL_c(jj)=(d_c(jj-1)-d_c(jj))/Mq_c(jj-1);
    %end
    
end

figure
plot(KL)
figure
plot(KL_c)