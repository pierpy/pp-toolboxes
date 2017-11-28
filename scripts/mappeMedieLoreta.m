clear all 
close all
clc

pathControlliLoretaF1 ='C:\Users\Pierpaolo\Documents\MATLAB\Microstati\Dati\MappeLoreta';

loretaControlliMapsF1=dir(pathControlliLoretaF1);

allMapsControlsF1=zeros(4,19);

for i= 3:size(loretaControlliMapsF1,1)
    x = load(strcat('C:\Users\Pierpaolo\Documents\MATLAB\Microstati\Dati\MappeLoreta\',loretaControlliMapsF1(i).name));
    
    allMapsControlsF1=[allMapsControlsF1;x];
   
end
allMapsControlsF1(1:4,:)=[];

figure
for k=1:size(x,1)
    subplot(2,2,k);
    topoplot(C(k,:), EEG.chanlocs)
    title('Loreta');
end