clear 
close all
clc

Nt=500;
Nm = 4;
Nch = 62;
% L = zeros(1,Nt);
% L(1:50)=1;
% L(51:100)=2;
% L(101:150)=3;
% L(151:end)=2;
V=zeros(Nch, Nt);
a = -1;
b = 1;
A = (b-a).*rand(Nm,Nt) + a;
T = (b-a).*rand(Nch,Nm) + a;
% [n_frame,n_chan] = size(T);
% h = eye(n_chan)-1/n_chan;
% T = T*h;	


for i=1:100
V(:,i)=T(:,1);
end
for i=101:200
V(:,i)=T(:,2);
end
for i=201:300
V(:,i)=T(:,3);
end
for i=301:500
V(:,i)=T(:,4);
end


data=V';

pmode=1;
average=0;
minclusters=1;
maxclusters=11;
ClusterNr = minclusters:maxclusters;
numclusters = length(ClusterNr);
[Ntr, Nchans]=size(data);


for q=1:numclusters
% Nmicrostates alg.
[maps,ind,assigned_cluster_ampl,exp_var] = NMicrostates(data,ClusterNr(q),pmode,20, average);
% calcolo CV
[NMicr.CV(q)] = CV(data,maps,ind,ClusterNr(q));
% calcolo W (dispersion measure)
[NMicr.W(q)] = W(data,ind,ClusterNr(q));
[NMicr.m(q)] = M(W(data,ind,ClusterNr(q)),Nchans,ClusterNr(q));
NMicr.template{q}=maps;
NMicr.clusters{q}=ind;
NMicr.amplitude{q}=assigned_cluster_ampl;
NMicr.exp_var(q)=exp_var;

%% simple k-means alg
% [clu, template]=kmeans(data,q);
% [Kmeans.CV(q)] = CV(data,template,clu,ClusterNr(q));
% [Kmeans.W(q)] = W(data,clu,ClusterNr(q));
% [Kmeans.m(q)] = M(Kmeans.W(q),Nchans,ClusterNr(q));
% 
% Kmeans.template{q}=template;
% Kmeans.clusters{q}=clu;

end


kl = KL(NMicr.m);
cv=NMicr.CV;
%kl = KL(Kmeans.m);
figure
subplot(211)
plot(kl);
title('KL')
subplot(212)
plot(cv);
title('CV')

Start=1;
Stop=size(data,1);

[MSClass,MSFit,gfp] = FitMicrostates(data',NMicr.template{1, 4},1,1,0.9,Start,Stop);
[OutFeature] = MSFeatures(MSClass,size(NMicr.template{1, 4},1),MSFit,gfp,1);

EEG.chanlocs = readlocs('coord.xyz');

figure;

for k=1:4
subplot(1,4,k);
topoplot(NMicr.template{1, 4}(k,:), EEG.chanlocs)
title(num2str(k))
end