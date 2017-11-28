clc
clear all
close all

%pathControlsMapsF1 ='C:\Users\Pierpaolo\Documents\MATLAB\Microstati\Dati\Results\SubjectFit\Mappe\Controlli\F1\4templates';
pathControlsMapsF2 ='C:\Users\Pierpaolo\Documents\MATLAB\Microstati\Dati\MappeMedie\groupSegmentation\Seg Controlf.KMeans\';
%pathLeftMapsF1 ='C:\Users\Pierpaolo\Documents\MATLAB\Microstati\Dati\Results\SubjectFit\Mappe\PazientiLeft\F1\4templates';
pathLeftMapsF2 ='C:\Users\Pierpaolo\Documents\MATLAB\Microstati\Dati\MappeMedie\groupSegmentation\Seg Leftf.KMeans\';
%pathRightMapsF1 ='C:\Users\Pierpaolo\Documents\MATLAB\Microstati\Dati\Results\SubjectFit\Mappe\PazientiRight\F1\4templates';
pathRightMapsF2 ='C:\Users\Pierpaolo\Documents\MATLAB\Microstati\Dati\MappeMedie\groupSegmentation\Seg Rightf.KMeans\';
pathPatient ='C:\Users\Pierpaolo\Documents\MATLAB\Microstati\Dati\MappeMedie\groupSegmentation\Seg Patient.KMeans\';

numeroTemplateControl='4';
numeroTemplateLeft='3';
numeroTemplateRight='4';
%numeroTemplatePatient='4';

%controlsMapsF1=dir(pathControlsMapsF1);
templateControl=strcat('Seg Controlf.KMeans.0',numeroTemplateControl,'.ep');
templateLeft=strcat('Seg Leftf.KMeans.0',numeroTemplateLeft,'.ep');
templateRight=strcat('Seg Rightf.KMeans.0',numeroTemplateRight,'.ep');
templatePatient=strcat('Seg Patient.KMeans.0',numeroTemplateRight,'.ep');




    %x = load(strcat('C:\Users\Pierpaolo\Documents\MATLAB\Microstati\Dati\Results\SubjectFit\Mappe\Controlli\F1\4templates\',controlsMapsF1(i).name));
    control = load(strcat(pathControlsMapsF2, templateControl));
    left = load(strcat(pathLeftMapsF2, templateLeft));
    right = load(strcat(pathRightMapsF2, templateRight));
%     patient = load(strcat(pathPatient, templatePatient));

%riordino e plot controlli
controlA=control(1,:);
controlB=control(4,:);
controlC=control(2,:);
controlD=control(3,:);
ordControl=[controlA;controlB;controlC;controlD];


EEG.chanlocs = readlocs('19ch2.xyz');



figure
for k=1:size(control,1)
    subplot(2,2,k);
    topoplot(control(k,:), EEG.chanlocs)
    title('Mean template Controls F2');
end
figure
for k=1:size(ordControl,1)
    subplot(2,2,k);
    topoplot(ordControl(k,:), EEG.chanlocs)
    title('Mean template Controls F2');
end

%riordino e plot left
leftA=left(2,:);
leftB=left(1,:);
leftD=left(4,:);
ordLeft=[leftA;leftB;leftD];
figure
for k=1:size(left,1)
    subplot(2,2,k);
    topoplot(left(k,:), EEG.chanlocs)
    title('Mean template PatientsLeft F2');
end
figure
for k=1:size(ordLeft,1)
    subplot(2,2,k);
    topoplot(ordLeft(k,:), EEG.chanlocs)
    title('Mean template PatientsLeft F2');
end

%riordino e plot right
rightA=right(1,:);
rightB=right(3,:);
rightC=right(2,:);
rightD=right(4,:);
ordRight=[rightA;rightB;rightC;rightD];
figure
for k=1:size(right,1)
    subplot(2,2,k);
    topoplot(right(k,:), EEG.chanlocs)
    title('Mean template PatientsRight F2');
end
figure
for k=1:size(ordRight,1)
    subplot(2,2,k);
    topoplot(ordRight(k,:), EEG.chanlocs)
    title('Mean template PatientsRight F2');
end

figure
for k=1:size(Result.Template{1, 3},2)
    subplot(2,2,k);
    topoplot(Result.Template{1, 3}(:,k), EEG.chanlocs)
    title('Mean template Patient F2');
end

figure
for k=1:size(Seg_GfpmaxAgagug_KMeans_04,1)
    subplot(2,2,k);
    topoplot(Seg_GfpmaxAgagug_KMeans_04(k,:), EEG.chanlocs)
    title('Mean template Patient F2');
end
saveeph('C:\Users\Pierpaolo\Documents\MATLAB\Microstati\Dati\MappeMedie\groupTemplatesControl.ep',ordControl)
saveeph('C:\Users\Pierpaolo\Documents\MATLAB\Microstati\Dati\MappeMedie\groupTemplatesLeft.ep',ordLeft)
saveeph('C:\Users\Pierpaolo\Documents\MATLAB\Microstati\Dati\MappeMedie\groupTemplatesRight.ep',ordRight)