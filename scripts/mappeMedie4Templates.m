clc
clear all
close all

path = 'C:\Users\Pierpaolo\Documents\MATLAB\MicrostatiBetaMotorio\rest_templates';
%controlsMapsF1=dir(pathControlsMapsF1);
templateFiles=dir(path);

allMaps=zeros(4,123);

for i= 3:size(templateFiles,1)
    %x = load(strcat('C:\Users\Pierpaolo\Documents\MATLAB\Microstati\Dati\Results\SubjectFit\Mappe\Controlli\F1\4templates\',controlsMapsF1(i).name));
    y = load(strcat('C:\Users\Pierpaolo\Documents\MATLAB\MicrostatiBetaMotorio\rest_templates\',templateFiles(i).name));
    %allMapsControlsF1=[allMapsControlsF1;x];
    allMaps=[allMaps;y];
end

allMaps(1:4,:)=[];
[IDX,templates] = kmeans(allMaps,4);
EEG.chanlocs = readlocs('coord.xyz');

figure
for k=1:size(templates,1)
    subplot(2,2,k);
    topoplot(templates(k,:), EEG.chanlocs)
    title('Mean template Rest');
end
