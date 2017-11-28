clear all
close all
clc

% pathGroup1=uigetdir;
% pathGroup2=uigetdir;
pathGroup1='C:\Users\Pierpaolo\Documents\MATLAB\Microstati\Dati\MappeIndividuali\healthy';
pathGroup2='C:\Users\Pierpaolo\Documents\MATLAB\Microstati\Dati\MappeIndividuali\right';
group1=dir(pathGroup1);group1(1)=[];group1(1)=[];
group2=dir(pathGroup2);group2(1)=[];group2(1)=[];

for k=1:size(group1) 
    G1=load(strcat(pathGroup1,'\',group1(k).name));
    AG1(k,:)=G1(1,:);
    BG1(k,:)=G1(2,:);
    CG1(k,:)=G1(3,:);
    DG1(k,:)=G1(4,:);
end
for j=1:size(group2)
    G2=load(strcat(pathGroup2,'\',group2(j).name));
    AG2(j,:)=G2(1,:);
    BG2(j,:)=G2(2,:);
    CG2(j,:)=G2(3,:);
    DG2(j,:)=G2(4,:);
end
% inverto le simmetriche a mano %%%%%
EEG.chanlocs = readlocs('19ch2.xyz');
for k=1:size(DG1,1)
%figure
subplot(3,10,k)
topoplot(DG1(k,:), EEG.chanlocs)
title(num2str(k))
end
DG1(2,:)=DG1(2,:).*-1;
DG1(8:10,:)=DG1(8:10,:).*-1;
DG1(6,:)=DG1(6,:).*-1;
DG1(16:17,:)=DG1(16:17,:).*-1;
DG1(14,:)=DG1(14,:).*-1;
DG2(22,:)=DG2(22,:).*-1;
DG2(27,:)=DG2(27,:).*-1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%