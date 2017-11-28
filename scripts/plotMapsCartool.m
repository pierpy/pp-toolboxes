clc
clear all
close all

%pathControlsMapsF1 ='C:\Users\Pierpaolo\Documents\MATLAB\Microstati\Dati\Results\SubjectFit\Mappe\Controlli\F1\4templates';
pathControlsMapsF2 ='C:\Users\Pierpaolo\Documents\MATLAB\Microstati\Dati\Results\SubjectFit\Mappe\Controlli\F2\4templates';
%pathLeftMapsF1 ='C:\Users\Pierpaolo\Documents\MATLAB\Microstati\Dati\Results\SubjectFit\Mappe\PazientiLeft\F1\4templates';
pathLeftMapsF2 ='C:\Users\Pierpaolo\Documents\MATLAB\Microstati\Dati\Results\SubjectFit\Mappe\PazientiLeft\F2\4templates';
%pathRightMapsF1 ='C:\Users\Pierpaolo\Documents\MATLAB\Microstati\Dati\Results\SubjectFit\Mappe\PazientiRight\F1\4templates';
pathRightMapsF2 ='C:\Users\Pierpaolo\Documents\MATLAB\Microstati\Dati\Results\SubjectFit\Mappe\PazientiRight\F2\4templates';


%controlsMapsF1=dir(pathControlsMapsF1);
controlsMapsF2=dir(pathControlsMapsF2);

%leftMapsF1=dir(pathLeftMapsF1);
leftMapsF2=dir(pathLeftMapsF2);

%rightMapsF1=dir(pathRightMapsF1);
rightMapsF2=dir(pathRightMapsF2);

%allMapsControlsF1=zeros(4,19);
allMapsControlsF2=zeros(4,19);
%allMapsLeftF1=zeros(4,19);
allMapsLeftF2=zeros(4,19);
%allMapsRightF1=zeros(4,19);
allMapsRightF2=zeros(4,19);

for i= 3:size(controlsMapsF2,1)
    %x = load(strcat('C:\Users\Pierpaolo\Documents\MATLAB\Microstati\Dati\Results\SubjectFit\Mappe\Controlli\F1\4templates\',controlsMapsF1(i).name));
    y = load(strcat('C:\Users\Pierpaolo\Documents\MATLAB\Microstati\Dati\Results\SubjectFit\Mappe\Controlli\F2\4templates\',controlsMapsF2(i).name));
    %allMapsControlsF1=[allMapsControlsF1;x];
    allMapsControlsF2=[allMapsControlsF2;y];
end

for j=3:size(leftMapsF2)
    %z = load(strcat('C:\Users\Pierpaolo\Documents\MATLAB\Microstati\Dati\Results\SubjectFit\Mappe\PazientiLeft\F1\4templates\',leftMapsF1(j).name));
    xx = load(strcat('C:\Users\Pierpaolo\Documents\MATLAB\Microstati\Dati\Results\SubjectFit\Mappe\PazientiLeft\F2\4templates\',leftMapsF2(j).name));
    %allMapsLeftF1=[allMapsLeftF1;z];
    allMapsLeftF2=[allMapsLeftF2;xx];
end

for k =3:size(rightMapsF2)
    %yy = load(strcat('C:\Users\Pierpaolo\Documents\MATLAB\Microstati\Dati\Results\SubjectFit\Mappe\PazientiRight\F1\4templates\',rightMapsF1(k).name));
    zz = load(strcat('C:\Users\Pierpaolo\Documents\MATLAB\Microstati\Dati\Results\SubjectFit\Mappe\PazientiRight\F2\4templates\',rightMapsF2(k).name));
    %allMapsRightF1=[allMapsRightF1;yy];
    allMapsRightF2=[allMapsRightF2;zz];
end

%allMapsControlsF1(1:4,:)=[];
allMapsControlsF2(1:4,:)=[];
%allMapsLeftF1(1:4,:)=[];
allMapsLeftF2(1:4,:)=[];
%allMapsRightF1(1:4,:)=[];
allMapsRightF2(1:4,:)=[];

%allMapsF1=[allMapsControlsF1;allMapsLeftF1;allMapsRightF1];
allMapsF2=[allMapsControlsF2;allMapsLeftF2;allMapsRightF2];
allPatient = [allMapsLeftF2;allMapsRightF2];
%[IDX_ControlsF1,templatesControlsF1] = kmeans(allMapsControlsF1,4);
[IDX_ControlsF2,templatesControlsF2] = kmeans(allMapsControlsF2,4);
%[IDX_LeftF1,templatesLeftF1] = kmeans(allMapsLeftF1,4);
[IDX_LeftF2,templatesLeftF2] = kmeans(allMapsLeftF2,4);
%[IDX_RightF1,templatesRightF1] = kmeans(allMapsRightF1,4);
[IDX_RightF2,templatesRightF2] = kmeans(allMapsRightF2,4);
%[IDX_F1,templatesF1] = kmeans(allMapsF1,4);
[IDX_F2,templatesF2] = kmeans(allMapsF2,4);



EEG.chanlocs = readlocs('19ch2.xyz');

% figure
% for k=1:size(templatesControlsF1,1)
%     subplot(2,2,k);
%     topoplot(templatesControlsF1(k,:), EEG.chanlocs)
%     title('Mean template Controls F1');
% end

figure
for k=1:size(templatesControlsF2,1)
    subplot(2,2,k);
    topoplot(templatesControlsF2(k,:), EEG.chanlocs)
    title('Mean template Controls F2');
end
% 
% figure
% for k=1:size(templatesControlsF1,1)
%     subplot(2,2,k);
%     topoplot(templatesLeftF1(k,:), EEG.chanlocs)
%     title('Mean template PatientsLeft F1');
% end

figure
for k=1:size(templatesControlsF2,1)
    subplot(2,2,k);
    topoplot(templatesLeftF2(k,:), EEG.chanlocs)
    title('Mean template PatientsLeft F2');
end
% 
% figure
% for k=1:size(templatesControlsF1,1)
%     subplot(2,2,k);
%     topoplot(templatesRightF1(k,:), EEG.chanlocs)
%     title('Mean template PatientsRight F1');
% end

figure
for k=1:size(templatesControlsF2,1)
    subplot(2,2,k);
    topoplot(templatesRightF2(k,:), EEG.chanlocs)
    title('Mean template PatientsRight F2');
end

% figure
% for k=1:size(templatesF1,1)
%     subplot(2,2,k);
%     topoplot(templatesF1(k,:), EEG.chanlocs)
%     title('Mean template F1');
% end

figure
for k=1:size(templatesF2,1)
    subplot(2,2,k);
    topoplot(templatesF2(k,:), EEG.chanlocs)
    title('Mean template F2');
end

% figure
% for k=1:size(xx,1)
%     subplot(2,2,k);
%     topoplot(xx(k,:), EEG.chanlocs)
%     title('Mean template F2');
% end
% figure
% for k=1:40
%     
%     subplot(10,4,k);
%     topoplot(allMapsControlsF2(k,:), EEG.chanlocs)
%     
%     
%     
%     
% %     title('tutte left');
% end
% figure
% for jj=k:size(allMapsControlsF2,1)
%      subplot(10,4,jj-39);
%     topoplot(allMapsControlsF2(jj,:), EEG.chanlocs)
% 
% end

