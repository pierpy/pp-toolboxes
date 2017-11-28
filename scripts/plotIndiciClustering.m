clear all
close all
clc

path='C:\Users\Pierpaolo\Documents\MATLAB\Microstati\Dati\Results\RESULT LAST\indiciKlustering\controlli\';

cd(path);
elencoFiles=dir;
elencoFiles(1)=[];elencoFiles(1)=[];

for k=1:size(elencoFiles,1)
filename = elencoFiles(k).name;
[cluster(:,k), maps(:,k), gev(:,k), cv(:,k), w(:,k), kl(:,k)] = textread(strcat(path, elencoFiles(k).name), '%f %f %f %f %f %f','headerlines', 2) ;
end
% plot(maps, kl)
% hold on
% plot()