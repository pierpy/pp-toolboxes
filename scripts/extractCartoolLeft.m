clear all
close all
clc


lettere={'A';'B';'C';'D';'E';'F';'G';'H';'I';'J';'K';'L';'M';'N';'O';'P';'Q';'R';'S';'T';'U';'V';'W';'Y';'Z';'AA';'AB';'AC';'AD';'AE';'AF'};  
labels={'meanGevA';'meanGevB';'meanGevC';'meanGevD';
        'meanOccA';'meanOccB';'meanOccC';'meanOccD';
        'meanTimeCovA';'meanTimeCovB';'meanTimeCovC';'meanTimeCovD';
        'meanNumTFA';'meanNnumTFB';'meanNnumTFC';'meanNnumTFD';
        'meanDurA';'meanDurB';'meanDurC';'meanDurD';
        'meanCorrA';'meanCorrB';'meanCorrC';'meanCorrD'};
    
filespath='C:\Users\Pierpaolo\Documents\MATLAB\Microstati\Dati\Results\groupWiseOrdinati\right\3mappe\xls_metriche';
cd(filespath);
elencoFiles=dir;
elencoFiles(1)=[];elencoFiles(1)=[];
%xlswrite('indici.xlsx',labels,'A1:A24')

for k=1:size(elencoFiles,1)
% filename=uigetfile('*.csv', 'Selezione il file csv');
filename=elencoFiles(k).name;
[num, txt, row]=xlsread(filename);

NumTF=num(:,1);
MeanDur=num(:,4)./1000000;

MeanCorr=row(:,5);
MeanCorr(1)=[];

for ii=1:size(MeanCorr,1)
    if ~iscellstr(MeanCorr(ii))
        %MeanCorr(ii)='0.0000';
        MeanCorr{ii}='0.000000';
        
    end
end

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
for ii=1:size(Occ,1)
    if ~iscellstr(Occ(ii))
        %MeanCorr(ii)='0.0000';
        Occ{ii}='0.000000';
        
    end
end
Occ=cell2mat(Occ);
Occ=str2num(Occ);


gev_a=Gev(1:3:end);gev_b=Gev(2:3:end);gev_d=Gev(3:3:end);

occ_a=Occ(1:3:end);occ_b=Occ(2:3:end);occ_d=Occ(3:3:end);

timeCov_a=TimeCov(1:3:end);timeCov_b=TimeCov(2:3:end);timeCov_d=TimeCov(3:3:end);

numTF_a=NumTF(1:3:end);numTF_b=NumTF(2:3:end);numTF_d=NumTF(3:3:end);
Dur_a=MeanDur(1:3:end);Dur_b=MeanDur(2:3:end);Dur_d=MeanDur(3:3:end);

Corr_a=MeanCorr(1:3:end);Corr_b=MeanCorr(2:3:end);Corr_d=MeanCorr(3:3:end);
%%
meanGevA=mean(gev_a);meanGevB=mean(gev_b);meanGevD=mean(gev_d);
meanOccA=mean(occ_a);meanOccB=mean(occ_b);meanOccD=mean(occ_d);
meanTimeCovA=mean(timeCov_a);meanTimeCovB=mean(timeCov_b);meanTimeCovD=mean(timeCov_d);
meanNumTFA=mean(numTF_a);meanNnumTFB=mean(numTF_b);meanNnumTFD=mean(numTF_d);
meanDurA=mean(Dur_a);meanDurB=mean(Dur_b);meanDurD=mean(Dur_d);
meanCorrA=mean(Corr_a);meanCorrB=mean(Corr_b);meanCorrD=mean(Corr_d);

risultati=[meanGevA;meanGevB;meanGevD;
           meanOccA;meanOccB;meanOccD;
           meanTimeCovA;meanTimeCovB;meanTimeCovD;
           meanNumTFA;meanNnumTFB;meanNnumTFD;
           meanDurA;meanDurB;meanDurD;
           meanCorrA;meanCorrB;meanCorrD];

        
    
culumn=strcat(lettere(k+1),'1:',lettere(k+1),'30');
culumn=cell2mat(culumn);
xlswrite('indici.xlsx',labels,'A1:A24')
xlswrite('indici.xlsx',risultati,culumn)
end

