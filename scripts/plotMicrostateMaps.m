clc
clear all
close all
prompt = 'sigla soggetto:';
dlg_title = 'Input';
answer = inputdlg(prompt,dlg_title);
filename = strcat(answer(1),'.mat');
pathfile = strcat('C:\Users\Pierpaolo\Documents\MATLAB\Microstati\Dati\Dati_occhi_aperti\segmentationResults\',strcat(answer(1),'.mat'));
load(pathfile{1});
microStateA = data.mappe(1,:);
microStateB = data.mappe(2,:);
microStateC = data.mappe(3,:);
microStateD = data.mappe(4,:);
EEG.chanlocs = readlocs('19ch2.xyz');
figure
topoplot(microStateA, EEG.chanlocs)
figure
topoplot(microStateB, EEG.chanlocs)
figure
topoplot(microStateC, EEG.chanlocs)
figure
topoplot(microStateD, EEG.chanlocs)