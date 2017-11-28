clear 
close all
clc
tic

folder_path = uigetdir('C:\Users\Pierpaolo\Documents\MATLAB\microstates_tool\dati');
S = {'ModKMeans', 'AAHC', 'SimpleKMeans'};
[s,v] = listdlg('PromptString','Select a file:',...
                'SelectionMode','single',...
                'ListString',S);
algorithm = S(s);
algorithm = algorithm{1, 1};
risposta = inputdlg({'High-pass filtering','Low-pass filtering','Frequencies to be excluded (in Hz)',...
                     'Decimation factor','Polarity Mode (0=no, 1=yes)',...
                     'Min Clusters (integer number)',...
                     'Max Clusters (integer number)',...
                     'Z-score (0 = no, 1 = yes)',...
                     'Subjects to analyze (integer)'},...
                     'Processing parameter selection',1,{'none','none','none','none','0','1','10','0','1'});

h_filtering = str2num(risposta{1});
l_filtering = str2num(risposta{2});
freq_notch = str2num(risposta{3});
dec=str2num(risposta{4});
pmode = str2num(risposta{5});
minclusters = str2num(risposta{6});
maxclusters = str2num(risposta{7});
zscore = str2num(risposta{8});
subjects_to_analyze = str2num(risposta{9});

fprintf('Riepilogo parametri:\n ')
fprintf(strcat('h_filtering: ', num2str(h_filtering),'\n'))
fprintf(strcat('l_filtering: ', num2str(l_filtering),'\n'))
fprintf(strcat('notch: ', num2str(freq_notch),'\n'))
fprintf(strcat('decimation factor: ', num2str(dec),'\n'))
fprintf(strcat('polarity mode: ', num2str(pmode),'\n'))
fprintf(strcat('min clusters: ', num2str(minclusters),'\n'))
fprintf(strcat('max clusters: ', num2str(maxclusters),'\n'))
fprintf(strcat('zscore: ', num2str(zscore),'\n'))
fprintf(strcat('Subjects to analyze: ', num2str(subjects_to_analyze),'\n'))


  if subjects_to_analyze == 'all'
            subjects_to_analyze = size(subjects,1);
  end
   for subject = 1:subjects_to_analyze     
    [SegmentationResult{subject}] = Microstates_Analysis(algorithm , subject, risposta, folder_path);
   end
toc

%% some inizializations
templates=[];
KL = [];
CV = [];
E_VAR = [];
n_clusters = 1;

EEG.chanlocs = readlocs('Layout_TMS.xyz');
 %final_templates=SegmentationResult{1, 1}.template{1, n_clusters};
% final_templates=SegmentationResult{1, 2}.template{1, 1};
%% evenutale orientamento delle mappe subject-wise  
% [orientend_maps]=orient_maps(SegmentationResult, nclusters, EEG);

%% calcolo templates di gruppo
for s = 1: size(SegmentationResult,2)
    templates = [templates; SegmentationResult{1,s}.template{1,n_clusters}];
end

figure;
for k=1:size(templates,1)
    subplot(9,8,k);
    topoplot(templates(k,:), EEG.chanlocs)
    title(num2str(k))
end

[final_templates] = mean_templates(templates, 4, 0);
n_clusters=4;
figure;
for k=1:n_clusters 
    subplot(1,n_clusters ,k);
    topoplot(final_templates(k,:), EEG.chanlocs)
    title(num2str(k))
end

%% backfitting group mean templates
for subject = 1:size(SegmentationResult,2)
    data = NMicr.original_allData;
    Start=1;
    Stop=size(data,1);  
    [MSClass,MSFit,gfp] = FitMicrostates(data',final_templates,1,1,0.9,Start,Stop);
    [OutFeature] = MSFeatures(MSClass,size(final_templates,1),MSFit,gfp,1);
    BackFitting{subject}=OutFeature;
    BackFitting{subject}.MSClass = MSClass;
    BackFitting{subject}.MSFit = MSFit;
    BackFitting{subject}.gfp = gfp;
    clearvars data OutFeature 
end

for i=1:size(SegmentationResult,2)
    durate(i,:)=BackFitting{1,i}.meandur;
end
% % 
 media = mean(durate)*(1000/125)


  

