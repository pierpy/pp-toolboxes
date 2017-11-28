clear all
close all
clc

    
filespath='C:\Users\Pierpaolo\Documents\MATLAB\Microstati\Dati\Results\groupWiseOrdinati\right\4mappe\xls_sintassi\20ripetizioni';
cd(filespath);
elencoFiles=dir;
elencoFiles(1)=[];elencoFiles(1)=[];
%xlswrite('indici.xlsx',labels,'A1:A24')

for k=1:size(elencoFiles,1)
% filename=uigetfile('*.csv', 'Selezione il file csv');
filename=elencoFiles(k).name;
[num, txt, row]=xlsread(filename);

for jj=2:17
prob=row(:,jj);
prob(1)=[];
for ii=1:size(prob,1)
    if ~iscellstr(prob(ii))
        %MeanCorr(ii)='0.0000';
        prob{ii}='0.000000';
        
    end
end
prob=cell2mat(prob);
prob=str2num(prob);
mean=mean(prob);
transizioni(k,jj-1)=mean;
clear prob mean
end



        
    
% culumn=strcat(lettere(k+1),'1:',lettere(k+1),'30');
% culumn=cell2mat(culumn);
% xlswrite('transizioni.xlsx',labels,'A1:A24')
%xlswrite('transizioni.xlsx',transizioni)
clear num txt row
end

