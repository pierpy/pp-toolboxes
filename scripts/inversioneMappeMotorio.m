clear all
close all
clc

pathGroup1='C:\Users\Pierpaolo\Documents\MATLAB\MicrostatiStimolazione\templates\Sham\Time-2\4templates';

group1=dir(pathGroup1);group1(1)=[];group1(1)=[];


for k=1:size(group1) 
    G1=load(strcat(pathGroup1,'\',group1(k).name));
    AG1(k,:)=G1(1,:);
    BG1(k,:)=G1(2,:);
    CG1(k,:)=G1(3,:);
    DG1(k,:)=G1(4,:);
%     EG1(k,:)=G1(5,:);
%      FG1(k,:)=G1(6,:);
end
% for j=1:size(group2)
%     G2=load(strcat(pathGroup2,'\',group2(j).name));
%     AG2(j,:)=G2(1,:);
%     BG2(j,:)=G2(2,:);
%     CG2(j,:)=G2(3,:);
%     DG2(j,:)=G2(4,:);
% end
% inverto le simmetriche a mano %%%%%
EEG.chanlocs = readlocs('coord.xyz');
figure;

for k=1:4
subplot(2,4,k);
topoplot(templates4(k,:), EEG.chanlocs)
title(num2str(k))
% title(num2str(group1(k).name(4:end-14)))
end
figure;

for k=1:size(BG1,1)
%figure
subplot(2,8,k)
topoplot(BG1(k,:), EEG.chanlocs)
title(num2str(k))
% title(num2str(group1(k).name(4:end-14)))
end
figure;

for k=1:size(CG1,1)
%figure
subplot(2,8,k)
topoplot(CG1(k,:), EEG.chanlocs)
title(num2str(k))
% title(num2str(group1(k).name(4:end-14)))
end
figure;

for k=1:size(DG1,1)
%figure
subplot(2,8,k)
topoplot(DG1(k,:), EEG.chanlocs)
title(num2str(k))
% title(num2str(group1(k).name(4:end-14)))
end
% figure;

% for k=1:size(EG1,1)
% subplot(2,8,k)
% topoplot(EG1(k,:), EEG.chanlocs)
% title(num2str(group1(k).name(4:end-14)))
% end

% figure;
% 
% for k=1:size(FG1,1)
% %figure
% subplot(3,6,k)
% topoplot(FG1(k,:), EEG.chanlocs)
% title(num2str(group1(k).name(4:end-14)))
% end



AG1(1,:)=AG1(1,:).*-1;
AG1(2,:)=AG1(2,:).*-1;
AG1(3,:)=AG1(3,:).*-1;
AG1(4,:)=AG1(4,:).*-1;
AG1(5,:)=AG1(5,:).*-1;
AG1(6,:)=AG1(6,:).*-1;
AG1(7,:)=AG1(7,:).*-1;
AG1(8,:)=AG1(8,:).*-1;
AG1(9,:)=AG1(9,:).*-1;
AG1(10,:)=AG1(10,:).*-1;
AG1(11,:)=AG1(11,:).*-1;
AG1(12,:)=AG1(12,:).*-1;
AG1(13,:)=AG1(13,:).*-1;
AG1(14,:)=AG1(14,:).*-1;
AG1(15,:)=AG1(15,:).*-1;
AG1(16,:)=AG1(16,:).*-1;
AG1(17,:)=AG1(17,:).*-1;
AG1(18,:)=AG1(18,:).*-1;

BG1(1,:)=BG1(1,:).*-1;
BG1(2,:)=BG1(2,:).*-1;
BG1(3,:)=BG1(3,:).*-1;
BG1(4,:)=BG1(4,:).*-1;
BG1(5,:)=BG1(5,:).*-1;
BG1(6,:)=BG1(6,:).*-1;
BG1(7,:)=BG1(7,:).*-1;
BG1(8,:)=BG1(8,:).*-1;
B1(9,:)= BG1(9,:).*-1;
BG1(10,:)=BG1(10,:).*-1;
BG1(11,:)=BG1(11,:).*-1;
BG1(12,:)=BG1(12,:).*-1;
BG1(13,:)=BG1(13,:).*-1;
BG1(14,:)=BG1(14,:).*-1;
BG1(15,:)=BG1(15,:).*-1;
BG1(16,:)=BG1(16,:).*-1;
BG1(17,:)=BG1(17,:).*-1;
BG1(18,:)=BG1(18,:).*-1;

CG1(1,:)=CG1(1,:).*-1;
CG1(2,:)=CG1(2,:).*-1;
CG1(3,:)=CG1(3,:).*-1;
CG1(4,:)=CG1(4,:).*-1;
CG1(5,:)=CG1(5,:).*-1;
CG1(6,:)=CG1(6,:).*-1;
CG1(7,:)=CG1(7,:).*-1;
CG1(8,:)=CG1(8,:).*-1;
CG1(9,:)= CG1(9,:).*-1;
CG1(10,:)=CG1(10,:).*-1;
CG1(11,:)=CG1(11,:).*-1;
CG1(12,:)=CG1(12,:).*-1;
CG1(13,:)=CG1(13,:).*-1;
CG1(14,:)=CG1(14,:).*-1;
CG1(15,:)=CG1(15,:).*-1;
CG1(16,:)=CG1(16,:).*-1;
CG1(17,:)=CG1(17,:).*-1;
CG1(18,:)=CG1(18,:).*-1;

DG1(1,:)=DG1(1,:).*-1;
DG1(2,:)=DG1(2,:).*-1;
DG1(3,:)=DG1(3,:).*-1;
DG1(4,:)=DG1(4,:).*-1;
DG1(5,:)=DG1(5,:).*-1;
DG1(6,:)=DG1(6,:).*-1;
DG1(7,:)=DG1(7,:).*-1;
DG1(8,:)=DG1(8,:).*-1;
DG1(9,:)= DG1(9,:).*-1;
DG1(10,:)=DG1(10,:).*-1;
DG1(11,:)=DG1(11,:).*-1;
DG1(12,:)=DG1(12,:).*-1;
DG1(13,:)=DG1(13,:).*-1;
DG1(14,:)=DG1(14,:).*-1;
DG1(15,:)=DG1(15,:).*-1;
DG1(16,:)=DG1(16,:).*-1;
DG1(17,:)=DG1(17,:).*-1;
DG1(18,:)=DG1(18,:).*-1;
% 
% % EG1(1,:)=EG1(1,:).*-1;
% % EG1(2,:)=EG1(2,:).*-1;
% % EG1(3,:)=EG1(3,:).*-1;
% % EG1(4,:)=EG1(4,:).*-1;
% % EG1(5,:)=EG1(5,:).*-1;
% % EG1(6,:)=EG1(6,:).*-1;
% % EG1(7,:)=EG1(7,:).*-1;
% % EG1(8,:)=EG1(8,:).*-1;
% % EG1(9,:)= EG1(9,:).*-1;
% % EG1(10,:)=EG1(10,:).*-1;
% % EG1(11,:)=EG1(11,:).*-1;
% % EG1(12,:)=EG1(12,:).*-1;
% % EG1(13,:)=EG1(13,:).*-1;
% % EG1(14,:)=EG1(14,:).*-1;
% % EG1(15,:)=EG1(15,:).*-1;
% % EG1(16,:)=EG1(16,:).*-1;
% % EG1(17,:)=EG1(17,:).*-1;
% % EG1(18,:)=EG1(18,:).*-1;
% % 
% % FG1(1,:)=FG1(1,:).*-1;
% % FG1(2,:)=FG1(2,:).*-1;
% % FG1(3,:)=FG1(3,:).*-1;
% % FG1(4,:)=FG1(4,:).*-1;
% % FG1(5,:)=FG1(5,:).*-1;
% % FG1(6,:)=FG1(6,:).*-1;
% % FG1(7,:)=FG1(7,:).*-1;
% % FG1(8,:)=FG1(8,:).*-1;
% % FG1(9,:)= FG1(9,:).*-1;
% % FG1(10,:)=FG1(10,:).*-1;
% % FG1(11,:)=FG1(11,:).*-1;
% % FG1(12,:)=FG1(12,:).*-1;
% % FG1(13,:)=FG1(13,:).*-1;
% % FG1(14,:)=FG1(14,:).*-1;
% % FG1(15,:)=FG1(15,:).*-1;
% % FG1(16,:)=FG1(16,:).*-1;
% % FG1(17,:)=FG1(17,:).*-1;
% % FG1(18,:)=FG1(18,:).*-1;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
allMaps = [AG1;BG1;CG1;DG1];

[IDX,templates] = kmeans(allMaps,4);

figure
for k=1:size(templates,1)
    subplot(2,4,k);
    topoplot(templates(k,:), EEG.chanlocs)
    title('Mean templates');
end
