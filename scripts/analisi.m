clear all
close all
clc


path='C:\Users\Pierpaolo\Documents\MATLAB\Microstati\Dati\concatenati\sani\';
files=dir(path);files(1)=[];files(1)=[];
settings.maxclusters=10;
settings.MaxIter=300;
% normalizzo %

for ii=1:size(files,2)
%data = load(strcat(path, files(ii).name));
data=load('C:\Users\Pierpaolo\Documents\MATLAB\Microstati\Dati\PazientiLeft\F2\AGAGUG_F2\All AgagugF.Export.txt');
data=data';
Result = lab_kmeans(data,settings);
risultati(ii)=Result;
end

load('C:\Users\Pierpaolo\Documents\MATLAB\Microstati\Dati\PazientiLeft\F2\AGAGUG_F2\Seg GfpmaxAgagug.KMeans\Seg GfpmaxAgagug.KMeans.04.ep');