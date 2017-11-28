clear all
close all
clc

labels={'meanGevA';'meanGevB';'meanGevC';'meanGevD';
        'meanOccA';'meanOccB';'meanOccC';'meanOccD';
        'meanTimeCovA';'meanTimeCovB';'meanTimeCovC';'meanTimeCovD';
        'meanNumTFA';'meanNnumTFB';'meanNnumTFC';'meanNnumTFD';
        'meanDurA';'meanDurB';'meanDurC';'meanDurD';
        'meanCorrA';'meanCorrB';'meanCorrC';'meanCorrD'};
    
pathControlliF1='C:\Users\Pierpaolo\Documents\MATLAB\Microstati\Dati\Results\SubjectFit\Metriche\F1\controlli';
pathLeftF1='C:\Users\Pierpaolo\Documents\MATLAB\Microstati\Dati\Results\SubjectFit\Metriche\F1\lesioniLeft';
pathRightF1='C:\Users\Pierpaolo\Documents\MATLAB\Microstati\Dati\Results\SubjectFit\Metriche\F1\lesioniRight';
pathControlliF2='C:\Users\Pierpaolo\Documents\MATLAB\Microstati\Dati\Results\SubjectFit\Metriche\F2\controlli';
pathLeftF2='C:\Users\Pierpaolo\Documents\MATLAB\Microstati\Dati\Results\SubjectFit\Metriche\F2\lesioniLeft';
pathRightF2='C:\Users\Pierpaolo\Documents\MATLAB\Microstati\Dati\Results\SubjectFit\Metriche\F2\lesioniRight';

%% Controlli F1
cd(pathControlliF1);
elencoFiles=dir;
elencoFiles(1)=[];elencoFiles(1)=[];


for k=1:size(elencoFiles,1)
% filename=uigetfile('*.csv', 'Selezione il file csv');
filename=elencoFiles(k).name;
[num, txt, row]=xlsread(filename);

NumTF=num(:,1);
MeanDur=num(:,4)./1000000;

MeanCorr=row(:,5);
MeanCorr(1)=[];
MeanCorr=cell2mat(MeanCorr);
MeanCorr=str2num(MeanCorr);

Gev=row(:,6);
Gev(1)=[];
Gev=cell2mat(Gev);
Gev=str2num(Gev);

TimeCov=row(:,8);
TimeCov(1)=[];
TimeCov=cell2mat(TimeCov);
TimeCov=str2num(TimeCov);

Occ=row(:,9);
Occ(1)=[];
Occ=cell2mat(Occ);
Occ=str2num(Occ);

gev_a=Gev(1:4:end);gev_b=Gev(2:4:end);gev_c=Gev(3:4:end);gev_d=Gev(4:4:end);

occ_a=Occ(1:4:end);occ_b=Occ(2:4:end);occ_c=Occ(3:4:end);occ_d=Occ(4:4:end);

timeCov_a=TimeCov(1:4:end);timeCov_b=TimeCov(2:4:end);timeCov_c=TimeCov(3:4:end);timeCov_d=TimeCov(4:4:end);

numTF_a=NumTF(1:4:end);numTF_b=NumTF(2:4:end);numTF_c=NumTF(3:4:end);numTF_d=NumTF(4:4:end);
Dur_a=MeanDur(1:4:end);Dur_b=MeanDur(2:4:end);Dur_c=MeanDur(3:4:end);Dur_d=MeanDur(4:4:end);

Corr_a=MeanCorr(1:4:end);Corr_b=MeanCorr(2:4:end);Corr_c=MeanCorr(3:4:end);Corr_d=MeanCorr(4:4:end);

meanGevA=mean(gev_a);meanGevB=mean(gev_b);meanGevC=mean(gev_c);meanGevD=mean(gev_d);
meanOccA=mean(occ_a);meanOccB=mean(occ_b);meanOccC=mean(occ_c);meanOccD=mean(occ_d);
meanTimeCovA=mean(timeCov_a);meanTimeCovB=mean(timeCov_b);meanTimeCovC=mean(timeCov_c);meanTimeCovD=mean(timeCov_d);
meanNumTFA=mean(numTF_a);meanNnumTFB=mean(numTF_b);meanNnumTFC=mean(numTF_c);meanNnumTFD=mean(numTF_d);
meanDurA=mean(Dur_a);meanDurB=mean(Dur_b);meanDurC=mean(Dur_c);meanDurD=mean(Dur_d);
meanCorrA=mean(Corr_a);meanCorrB=mean(Corr_b);meanCorrC=mean(Corr_c);meanCorrD=mean(Corr_d);

risultatiControlliF1(:,k)=[meanGevA;meanGevB;meanGevC;meanGevD;
           meanOccA;meanOccB;meanOccC;meanOccD;
           meanTimeCovA;meanTimeCovB;meanTimeCovC;meanTimeCovD;
           meanNumTFA;meanNnumTFB;meanNnumTFC;meanNnumTFD;
           meanDurA;meanDurB;meanDurC;meanDurD;
           meanCorrA;meanCorrB;meanCorrC;meanCorrD];
end
clear elencoFiles;

%% Left F1
cd(pathLeftF1);
elencoFiles=dir;
elencoFiles(1)=[];elencoFiles(1)=[];


for k=1:size(elencoFiles,1)
% filename=uigetfile('*.csv', 'Selezione il file csv');
filename=elencoFiles(k).name;
[num, txt, row]=xlsread(filename);

NumTF=num(:,1);
MeanDur=num(:,4)./1000000;

MeanCorr=row(:,5);
MeanCorr(1)=[];
MeanCorr=cell2mat(MeanCorr);
MeanCorr=str2num(MeanCorr);

Gev=row(:,6);
Gev(1)=[];
Gev=cell2mat(Gev);
Gev=str2num(Gev);

TimeCov=row(:,8);
TimeCov(1)=[];
TimeCov=cell2mat(TimeCov);
TimeCov=str2num(TimeCov);

Occ=row(:,9);
Occ(1)=[];
Occ=cell2mat(Occ);
Occ=str2num(Occ);

gev_a=Gev(1:4:end);gev_b=Gev(2:4:end);gev_c=Gev(3:4:end);gev_d=Gev(4:4:end);

occ_a=Occ(1:4:end);occ_b=Occ(2:4:end);occ_c=Occ(3:4:end);occ_d=Occ(4:4:end);

timeCov_a=TimeCov(1:4:end);timeCov_b=TimeCov(2:4:end);timeCov_c=TimeCov(3:4:end);timeCov_d=TimeCov(4:4:end);

numTF_a=NumTF(1:4:end);numTF_b=NumTF(2:4:end);numTF_c=NumTF(3:4:end);numTF_d=NumTF(4:4:end);
Dur_a=MeanDur(1:4:end);Dur_b=MeanDur(2:4:end);Dur_c=MeanDur(3:4:end);Dur_d=MeanDur(4:4:end);

Corr_a=MeanCorr(1:4:end);Corr_b=MeanCorr(2:4:end);Corr_c=MeanCorr(3:4:end);Corr_d=MeanCorr(4:4:end);

meanGevA=mean(gev_a);meanGevB=mean(gev_b);meanGevC=mean(gev_c);meanGevD=mean(gev_d);
meanOccA=mean(occ_a);meanOccB=mean(occ_b);meanOccC=mean(occ_c);meanOccD=mean(occ_d);
meanTimeCovA=mean(timeCov_a);meanTimeCovB=mean(timeCov_b);meanTimeCovC=mean(timeCov_c);meanTimeCovD=mean(timeCov_d);
meanNumTFA=mean(numTF_a);meanNnumTFB=mean(numTF_b);meanNnumTFC=mean(numTF_c);meanNnumTFD=mean(numTF_d);
meanDurA=mean(Dur_a);meanDurB=mean(Dur_b);meanDurC=mean(Dur_c);meanDurD=mean(Dur_d);
meanCorrA=mean(Corr_a);meanCorrB=mean(Corr_b);meanCorrC=mean(Corr_c);meanCorrD=mean(Corr_d);

risultatiLeftF1(:,k)=[meanGevA;meanGevB;meanGevC;meanGevD;
           meanOccA;meanOccB;meanOccC;meanOccD;
           meanTimeCovA;meanTimeCovB;meanTimeCovC;meanTimeCovD;
           meanNumTFA;meanNnumTFB;meanNnumTFC;meanNnumTFD;
           meanDurA;meanDurB;meanDurC;meanDurD;
           meanCorrA;meanCorrB;meanCorrC;meanCorrD];
end
clear elencoFiles;
%% Right F1 
cd(pathRightF1);
elencoFiles=dir;
elencoFiles(1)=[];elencoFiles(1)=[];


for k=1:size(elencoFiles,1)
% filename=uigetfile('*.csv', 'Selezione il file csv');
filename=elencoFiles(k).name;
[num, txt, row]=xlsread(filename);

NumTF=num(:,1);
MeanDur=num(:,4)./1000000;

MeanCorr=row(:,5);
MeanCorr(1)=[];
MeanCorr=cell2mat(MeanCorr);
MeanCorr=str2num(MeanCorr);

Gev=row(:,6);
Gev(1)=[];
Gev=cell2mat(Gev);
Gev=str2num(Gev);

TimeCov=row(:,8);
TimeCov(1)=[];
TimeCov=cell2mat(TimeCov);
TimeCov=str2num(TimeCov);

Occ=row(:,9);
Occ(1)=[];
Occ=cell2mat(Occ);
Occ=str2num(Occ);

gev_a=Gev(1:4:end);gev_b=Gev(2:4:end);gev_c=Gev(3:4:end);gev_d=Gev(4:4:end);

occ_a=Occ(1:4:end);occ_b=Occ(2:4:end);occ_c=Occ(3:4:end);occ_d=Occ(4:4:end);

timeCov_a=TimeCov(1:4:end);timeCov_b=TimeCov(2:4:end);timeCov_c=TimeCov(3:4:end);timeCov_d=TimeCov(4:4:end);

numTF_a=NumTF(1:4:end);numTF_b=NumTF(2:4:end);numTF_c=NumTF(3:4:end);numTF_d=NumTF(4:4:end);
Dur_a=MeanDur(1:4:end);Dur_b=MeanDur(2:4:end);Dur_c=MeanDur(3:4:end);Dur_d=MeanDur(4:4:end);

Corr_a=MeanCorr(1:4:end);Corr_b=MeanCorr(2:4:end);Corr_c=MeanCorr(3:4:end);Corr_d=MeanCorr(4:4:end);

meanGevA=mean(gev_a);meanGevB=mean(gev_b);meanGevC=mean(gev_c);meanGevD=mean(gev_d);
meanOccA=mean(occ_a);meanOccB=mean(occ_b);meanOccC=mean(occ_c);meanOccD=mean(occ_d);
meanTimeCovA=mean(timeCov_a);meanTimeCovB=mean(timeCov_b);meanTimeCovC=mean(timeCov_c);meanTimeCovD=mean(timeCov_d);
meanNumTFA=mean(numTF_a);meanNnumTFB=mean(numTF_b);meanNnumTFC=mean(numTF_c);meanNnumTFD=mean(numTF_d);
meanDurA=mean(Dur_a);meanDurB=mean(Dur_b);meanDurC=mean(Dur_c);meanDurD=mean(Dur_d);
meanCorrA=mean(Corr_a);meanCorrB=mean(Corr_b);meanCorrC=mean(Corr_c);meanCorrD=mean(Corr_d);

risultatiRightF1(:,k)=[meanGevA;meanGevB;meanGevC;meanGevD;
           meanOccA;meanOccB;meanOccC;meanOccD;
           meanTimeCovA;meanTimeCovB;meanTimeCovC;meanTimeCovD;
           meanNumTFA;meanNnumTFB;meanNnumTFC;meanNnumTFD;
           meanDurA;meanDurB;meanDurC;meanDurD;
           meanCorrA;meanCorrB;meanCorrC;meanCorrD];
end
clear elencoFiles;

%% Controlli F2
cd(pathControlliF2);
elencoFiles=dir;
elencoFiles(1)=[];elencoFiles(1)=[];


for k=1:size(elencoFiles,1)
% filename=uigetfile('*.csv', 'Selezione il file csv');
filename=elencoFiles(k).name;
[num, txt, row]=xlsread(filename);

NumTF=num(:,1);
MeanDur=num(:,4)./1000000;

MeanCorr=row(:,5);
MeanCorr(1)=[];
MeanCorr=cell2mat(MeanCorr);
MeanCorr=str2num(MeanCorr);

Gev=row(:,6);
Gev(1)=[];
Gev=cell2mat(Gev);
Gev=str2num(Gev);

TimeCov=row(:,8);
TimeCov(1)=[];
TimeCov=cell2mat(TimeCov);
TimeCov=str2num(TimeCov);

Occ=row(:,9);
Occ(1)=[];
Occ=cell2mat(Occ);
Occ=str2num(Occ);

gev_a=Gev(1:4:end);gev_b=Gev(2:4:end);gev_c=Gev(3:4:end);gev_d=Gev(4:4:end);

occ_a=Occ(1:4:end);occ_b=Occ(2:4:end);occ_c=Occ(3:4:end);occ_d=Occ(4:4:end);

timeCov_a=TimeCov(1:4:end);timeCov_b=TimeCov(2:4:end);timeCov_c=TimeCov(3:4:end);timeCov_d=TimeCov(4:4:end);

numTF_a=NumTF(1:4:end);numTF_b=NumTF(2:4:end);numTF_c=NumTF(3:4:end);numTF_d=NumTF(4:4:end);
Dur_a=MeanDur(1:4:end);Dur_b=MeanDur(2:4:end);Dur_c=MeanDur(3:4:end);Dur_d=MeanDur(4:4:end);

Corr_a=MeanCorr(1:4:end);Corr_b=MeanCorr(2:4:end);Corr_c=MeanCorr(3:4:end);Corr_d=MeanCorr(4:4:end);

meanGevA=mean(gev_a);meanGevB=mean(gev_b);meanGevC=mean(gev_c);meanGevD=mean(gev_d);
meanOccA=mean(occ_a);meanOccB=mean(occ_b);meanOccC=mean(occ_c);meanOccD=mean(occ_d);
meanTimeCovA=mean(timeCov_a);meanTimeCovB=mean(timeCov_b);meanTimeCovC=mean(timeCov_c);meanTimeCovD=mean(timeCov_d);
meanNumTFA=mean(numTF_a);meanNnumTFB=mean(numTF_b);meanNnumTFC=mean(numTF_c);meanNnumTFD=mean(numTF_d);
meanDurA=mean(Dur_a);meanDurB=mean(Dur_b);meanDurC=mean(Dur_c);meanDurD=mean(Dur_d);
meanCorrA=mean(Corr_a);meanCorrB=mean(Corr_b);meanCorrC=mean(Corr_c);meanCorrD=mean(Corr_d);

risultatiControlsF2(:,k)=[meanGevA;meanGevB;meanGevC;meanGevD;
           meanOccA;meanOccB;meanOccC;meanOccD;
           meanTimeCovA;meanTimeCovB;meanTimeCovC;meanTimeCovD;
           meanNumTFA;meanNnumTFB;meanNnumTFC;meanNnumTFD;
           meanDurA;meanDurB;meanDurC;meanDurD;
           meanCorrA;meanCorrB;meanCorrC;meanCorrD];
end
clear elencoFiles;
%% Left F2
cd(pathLeftF2);
elencoFiles=dir;
elencoFiles(1)=[];elencoFiles(1)=[];


for k=1:size(elencoFiles,1)
% filename=uigetfile('*.csv', 'Selezione il file csv');
filename=elencoFiles(k).name;
[num, txt, row]=xlsread(filename);

NumTF=num(:,1);
MeanDur=num(:,4)./1000000;

MeanCorr=row(:,5);
MeanCorr(1)=[];
MeanCorr=cell2mat(MeanCorr);
MeanCorr=str2num(MeanCorr);

Gev=row(:,6);
Gev(1)=[];
Gev=cell2mat(Gev);
Gev=str2num(Gev);

TimeCov=row(:,8);
TimeCov(1)=[];
TimeCov=cell2mat(TimeCov);
TimeCov=str2num(TimeCov);

Occ=row(:,9);
Occ(1)=[];
Occ=cell2mat(Occ);
Occ=str2num(Occ);

gev_a=Gev(1:4:end);gev_b=Gev(2:4:end);gev_c=Gev(3:4:end);gev_d=Gev(4:4:end);

occ_a=Occ(1:4:end);occ_b=Occ(2:4:end);occ_c=Occ(3:4:end);occ_d=Occ(4:4:end);

timeCov_a=TimeCov(1:4:end);timeCov_b=TimeCov(2:4:end);timeCov_c=TimeCov(3:4:end);timeCov_d=TimeCov(4:4:end);

numTF_a=NumTF(1:4:end);numTF_b=NumTF(2:4:end);numTF_c=NumTF(3:4:end);numTF_d=NumTF(4:4:end);
Dur_a=MeanDur(1:4:end);Dur_b=MeanDur(2:4:end);Dur_c=MeanDur(3:4:end);Dur_d=MeanDur(4:4:end);

Corr_a=MeanCorr(1:4:end);Corr_b=MeanCorr(2:4:end);Corr_c=MeanCorr(3:4:end);Corr_d=MeanCorr(4:4:end);

meanGevA=mean(gev_a);meanGevB=mean(gev_b);meanGevC=mean(gev_c);meanGevD=mean(gev_d);
meanOccA=mean(occ_a);meanOccB=mean(occ_b);meanOccC=mean(occ_c);meanOccD=mean(occ_d);
meanTimeCovA=mean(timeCov_a);meanTimeCovB=mean(timeCov_b);meanTimeCovC=mean(timeCov_c);meanTimeCovD=mean(timeCov_d);
meanNumTFA=mean(numTF_a);meanNnumTFB=mean(numTF_b);meanNnumTFC=mean(numTF_c);meanNnumTFD=mean(numTF_d);
meanDurA=mean(Dur_a);meanDurB=mean(Dur_b);meanDurC=mean(Dur_c);meanDurD=mean(Dur_d);
meanCorrA=mean(Corr_a);meanCorrB=mean(Corr_b);meanCorrC=mean(Corr_c);meanCorrD=mean(Corr_d);

risultatiLeftF2(:,k)=[meanGevA;meanGevB;meanGevC;meanGevD;
           meanOccA;meanOccB;meanOccC;meanOccD;
           meanTimeCovA;meanTimeCovB;meanTimeCovC;meanTimeCovD;
           meanNumTFA;meanNnumTFB;meanNnumTFC;meanNnumTFD;
           meanDurA;meanDurB;meanDurC;meanDurD;
           meanCorrA;meanCorrB;meanCorrC;meanCorrD];
end
clear elencoFiles;

%% Right F2
cd(pathRightF2);
elencoFiles=dir;
elencoFiles(1)=[];elencoFiles(1)=[];


for k=1:size(elencoFiles,1)
% filename=uigetfile('*.csv', 'Selezione il file csv');
filename=elencoFiles(k).name;
[num, txt, row]=xlsread(filename);

NumTF=num(:,1);
MeanDur=num(:,4)./1000000;

MeanCorr=row(:,5);
MeanCorr(1)=[];
MeanCorr=cell2mat(MeanCorr);
MeanCorr=str2num(MeanCorr);

Gev=row(:,6);
Gev(1)=[];
Gev=cell2mat(Gev);
Gev=str2num(Gev);

TimeCov=row(:,8);
TimeCov(1)=[];
TimeCov=cell2mat(TimeCov);
TimeCov=str2num(TimeCov);

Occ=row(:,9);
Occ(1)=[];
Occ=cell2mat(Occ);
Occ=str2num(Occ);

gev_a=Gev(1:4:end);gev_b=Gev(2:4:end);gev_c=Gev(3:4:end);gev_d=Gev(4:4:end);

occ_a=Occ(1:4:end);occ_b=Occ(2:4:end);occ_c=Occ(3:4:end);occ_d=Occ(4:4:end);

timeCov_a=TimeCov(1:4:end);timeCov_b=TimeCov(2:4:end);timeCov_c=TimeCov(3:4:end);timeCov_d=TimeCov(4:4:end);

numTF_a=NumTF(1:4:end);numTF_b=NumTF(2:4:end);numTF_c=NumTF(3:4:end);numTF_d=NumTF(4:4:end);
Dur_a=MeanDur(1:4:end);Dur_b=MeanDur(2:4:end);Dur_c=MeanDur(3:4:end);Dur_d=MeanDur(4:4:end);

Corr_a=MeanCorr(1:4:end);Corr_b=MeanCorr(2:4:end);Corr_c=MeanCorr(3:4:end);Corr_d=MeanCorr(4:4:end);

meanGevA=mean(gev_a);meanGevB=mean(gev_b);meanGevC=mean(gev_c);meanGevD=mean(gev_d);
meanOccA=mean(occ_a);meanOccB=mean(occ_b);meanOccC=mean(occ_c);meanOccD=mean(occ_d);
meanTimeCovA=mean(timeCov_a);meanTimeCovB=mean(timeCov_b);meanTimeCovC=mean(timeCov_c);meanTimeCovD=mean(timeCov_d);
meanNumTFA=mean(numTF_a);meanNnumTFB=mean(numTF_b);meanNnumTFC=mean(numTF_c);meanNnumTFD=mean(numTF_d);
meanDurA=mean(Dur_a);meanDurB=mean(Dur_b);meanDurC=mean(Dur_c);meanDurD=mean(Dur_d);
meanCorrA=mean(Corr_a);meanCorrB=mean(Corr_b);meanCorrC=mean(Corr_c);meanCorrD=mean(Corr_d);

risultatiRightF2(:,k)=[meanGevA;meanGevB;meanGevC;meanGevD;
           meanOccA;meanOccB;meanOccC;meanOccD;
           meanTimeCovA;meanTimeCovB;meanTimeCovC;meanTimeCovD;
           meanNumTFA;meanNnumTFB;meanNnumTFC;meanNnumTFD;
           meanDurA;meanDurB;meanDurC;meanDurD;
           meanCorrA;meanCorrB;meanCorrC;meanCorrD];
end
%% ttests
% for k=1:size(risultatiControlliF1,1)
%     [h_Controls_vs_LeftF1(k),p_Control_vs_LeftF1(k)]=ttest2(risultatiControlliF1(k,:),risultatiLeftF1(k,:));
% end
% for k=1:size(risultatiControlliF1,1)
%     [h_Controls_vs_RightF1(k),p_Control_vs_RightF1(k)]=ttest2(risultatiControlliF1(k,:),risultatiRightF1(k,:));
% end
% for k=1:size(risultatiControlliF1,1)
%     [h_Left_vs_RightF1(k),p_Left_vs_RightF1(k)]=ttest2(risultatiLeftF1(k,:),risultatiRightF1(k,:));
% end
% 
% for k=1:size(risultatiControlliF1,1)
%     [h_Controls_vs_LeftF2(k),p_Control_vs_LeftF2(k)]=ttest2(risultatiControlsF2(k,:),risultatiLeftF2(k,:));
% end
% for k=1:size(risultatiControlliF1,1)
%     [h_Controls_vs_RightF2(k),p_Control_vs_RightF2(k)]=ttest2(risultatiControlsF2(k,:),risultatiRightF2(k,:));
% end
% for k=1:size(risultatiControlliF1,1)
%     [h_Left_vs_RightF2(k),p_Left_vs_RightF2(k)]=ttest2(risultatiLeftF2(k,:),risultatiRightF2(k,:));
% end
