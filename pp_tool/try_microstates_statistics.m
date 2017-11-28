clear 
close all
clc
%files_path='/Users/pp/Documents/MATLAB/microstates/dati/sidfabmax';
folder_path='C:\Users\Pierpaolo\Documents\MATLAB\microstates_tool\dati\mov1';
% load(path);
cd(folder_path);
subjects = dir;
subjects(1)=[];
subjects(1)=[];
tic
for subject = 1:size(subjects,1)
    path = strcat(folder_path,'\', subjects(subject).name);
    cd(path)
    files=dir('*.interpolate.ep');
    total_data=[];

    for k=1:size(files,1)
        segment = load(strcat(path,'/',files(k).name));
        total_data=[total_data;segment];   
    end
   
    % any filtering
%     [total_data] = eegfilt(total_data',500,8,12);
%     total_data=total_data';
     % z score
    total_data = zscore(total_data);
    %  downsample
     total_data = downsample(total_data,5);
    % GFP maxima
    data = globalFieldPower(total_data');

    data=data';

    pmode=1;
    minclusters=1;
    maxclusters=12;
    ClusterNr = minclusters:maxclusters;
    numclusters = length(ClusterNr);
    [Ntr, Nchans]=size(data);


    for q=1:numclusters
    %% Nmicrostates alg.
    [maps,ind,assigned_cluster_ampl,exp_var] = NMicrostates(data,ClusterNr(q),pmode,20);
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
    NMicr.kl  =  KL(NMicr.m);
    total{subject}=NMicr;
    clearvars files segment
    cd(folder_path);
end
toc
% kl = KL(NMicr.m);
% 
% cv=NMicr.CV;
% %kl = KL(Kmeans.m);
% figure
% subplot(211)
% plot(kl);
% title('KL')
% subplot(212)
% plot(cv);
% title('CV')
% 

% 
% 
% 
EEG.chanlocs = readlocs('coord.xyz');


template=[];
for i=1:size(subjects,1)
    %KL(i,:)=total{i}.kl;
   template=[template;total{1, i}.template{1, 4}];
end

%orienta templates
inverted = [];
for  k=1:size(template,1)
      G1=total{1, k}.template{1, 4};
        
        figure;
        for kk=1:size(G1,1)
        subplot(1,4,kk);
        topoplot(G1(kk,:), EEG.chanlocs)
        title(num2str(kk))
        end
        
        button1 = questdlg('invertire mappa 1?');
        if strcmp(button1, 'Yes')
            disp('inverto mappa1');
            G1(1,:)=G1(1,:).*-1;
        end
        button2 = questdlg('invertire mappa 2?');
        if strcmp(button2, 'Yes')
            disp('inverto mappa2');
             G1(2,:)=G1(2,:).*-1;
        end
        button3 = questdlg('invertire mappa 3?');
        if strcmp(button3, 'Yes')
            disp('inverto mappa3');
            G1(3,:)=G1(3,:).*-1;
        end
        button4 = questdlg('invertire mappa 4?');
        if strcmp(button4, 'Yes')
            disp('inverto mappa4');
            G1(4,:)=G1(4,:).*-1;
        end
         
        close all
        figure;
        for kk=1:size(G1,1)
        subplot(1,4,kk);
        topoplot(G1(kk,:), EEG.chanlocs)
        title(num2str(kk))
        end
        
      inverted = [inverted;G1];
      
        clearvars G1
end

Nm = 4;
[final_templates, stab, V]=mean_templates(inverted, Nm);

figure;

for k=1:4
subplot(1,4,k);
topoplot(final_templates(k,:), EEG.chanlocs)
title(num2str(k))
end

%backfitting group mean templates
for subject = 1:size(subjects,1)
    path = strcat(folder_path,'\', subjects(subject).name);
    cd(path)
    files=dir('*.interpolate.ep');
    data=[];
    
    for kk=1:size(files,1)
        segment = load(strcat(path,'/',files(kk).name));
        data=[data;segment]; 
    end
    
    % any filtering
%      [data] = eegfilt(data',500,8,12);
%      data=data';
     data = zscore(data);
    %  downsample
    data_dnw = downsample(data,5);
    data_ok = globalFieldPower(data_dnw');
    Start=1;
    Stop=size(data_ok,2);

    [MSClass,MSFit,gfp] = FitMicrostates(data_ok,final_templates,1,1,0.9,Start,Stop);
    [OutFeature] = MSFeatures(MSClass,size(final_templates,1),MSFit,gfp,1);
    out{subject}=OutFeature;
    clearvars data OutFeature data_ok data_dnw
end


for i=1:11
    durate(i,:)=out{1,i}.meandur;
end

media = mean(durate)*10


  

