clear all
close all
clc

labels={'meanGevA';'meanGevB';'meanGevC';'meanGevD';
        'meanOccA';'meanOccB';'meanOccC';'meanOccD';
        'meanTimeCovA';'meanTimeCovB';'meanTimeCovC';'meanTimeCovD';
        'meanNumTFA';'meanNnumTFB';'meanNnumTFC';'meanNnumTFD';
        'meanDurA';'meanDurB';'meanDurC';'meanDurD';
        'meanCorrA';'meanCorrB';'meanCorrC';'meanCorrD'};
    
pathControls='C:\Users\Pierpaolo\Documents\MATLAB\Microstati\Dati\Results\GlobalFitCotrolsF2\csvFiles';
pathLeft='C:\Users\Pierpaolo\Documents\MATLAB\Microstati\Dati\Results\GlobalFitLeftF2\csvFiles';
pathRight='C:\Users\Pierpaolo\Documents\MATLAB\Microstati\Dati\Results\GlobalFitRightF2\csvFiles';

%% Global Controls F2
cd(pathControls);
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

globalControlsF2(:,k)=[meanGevA;meanGevB;meanGevC;meanGevD;
           meanOccA;meanOccB;meanOccC;meanOccD;
           meanTimeCovA;meanTimeCovB;meanTimeCovC;meanTimeCovD;
           meanNumTFA;meanNnumTFB;meanNnumTFC;meanNnumTFD;
           meanDurA;meanDurB;meanDurC;meanDurD;
           meanCorrA;meanCorrB;meanCorrC;meanCorrD];
end
clear elencoFiles;

%% Global Left F2
cd(pathLeft);
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

globalLeftF2(:,k)=[meanGevA;meanGevB;meanGevC;meanGevD;
           meanOccA;meanOccB;meanOccC;meanOccD;
           meanTimeCovA;meanTimeCovB;meanTimeCovC;meanTimeCovD;
           meanNumTFA;meanNnumTFB;meanNnumTFC;meanNnumTFD;
           meanDurA;meanDurB;meanDurC;meanDurD;
           meanCorrA;meanCorrB;meanCorrC;meanCorrD];
end
clear elencoFiles;
%% Global Right F2 
cd(pathRight);
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

globalRightF2(:,k)=[meanGevA;meanGevB;meanGevC;meanGevD;
           meanOccA;meanOccB;meanOccC;meanOccD;
           meanTimeCovA;meanTimeCovB;meanTimeCovC;meanTimeCovD;
           meanNumTFA;meanNnumTFB;meanNnumTFC;meanNnumTFD;
           meanDurA;meanDurB;meanDurC;meanDurD;
           meanCorrA;meanCorrB;meanCorrC;meanCorrD];
end



for k=1:size(globalControlsF2,1)
    [h_Controls_vs_Left(k),p_Control_vs_Left(k)]=ttest2(globalControlsF2(k,:),globalLeftF2(k,:));
end
for k=1:size(globalControlsF2,1)
    [h_Controls_vs_Right(k),p_Control_vs_Right(k)]=ttest2(globalControlsF2(k,:),globalRightF2(k,:));
end
for k=1:size(globalControlsF2,1)
    [h_Left_vs_Right(k),p_Left_vs_Right(k)]=ttest2(globalLeftF2(k,:),globalRightF2(k,:));
end