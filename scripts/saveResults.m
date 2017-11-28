clc
clear all
close all
prompt = 'sigla soggetto:';
dlg_title = 'Input';

answer = inputdlg(prompt,dlg_title);
filename=strcat(answer(1),'-04-mStParams(#FilesX#Params).txt');
%filename = uigetfile('*.txt', 'carica il file con i paramentri risulatanti analis');
cd('C:\Users\Pierpaolo\Documents\MATLAB\Microstati\Dati\Dati_occhi_aperti\segmentationResults');

data = importdata(filename{1});
params=data.data;
%labelsname=uigetfile('*.txt', 'carica i labels');
labelsname=strcat(answer(1),'-04-Labels(#TFx#Files).txt');
labels=importdata(labelsname{1});
%probname=uigetfile('*.txt', 'carica le probailità di transizione');
probname=strcat(answer(1),'-04-MicrostateProbabilities(#TFx#MicroSts).txt');
probabilities=importdata(probname{1});
%mappename=uigetfile('*.txt', 'carica le mappe');
mappename=strcat(answer(1),'-4-uStMaps-Temp.txt');
mappe=importdata(mappename{1});
data.params=params;
data.labels=labels;
data.prob=probabilities;
data.mappe=mappe;
save(answer{1},'data')




