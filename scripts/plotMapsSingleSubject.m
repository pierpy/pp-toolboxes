clc
clear all
close all

%pathControlsMapsF1 ='C:\Users\Pierpaolo\Documents\MATLAB\Microstati\Dati\Results\SubjectFit\Mappe\Controlli\F1\4templates';
singleMapControl ='C:\Users\Pierpaolo\Documents\MATLAB\Microstati\Dati\Results\SubjectFit\Mappe\Controlli\F2\4templates\';
%pathLeftMapsF1 ='C:\Users\Pierpaolo\Documents\MATLAB\Microstati\Dati\Results\SubjectFit\Mappe\PazientiLeft\F1\4templates';
singleMapLeft ='C:\Users\Pierpaolo\Documents\MATLAB\Microstati\Dati\Results\SubjectFit\Mappe\PazientiLeft\F2\4templates';
%pathRightMapsF1 ='C:\Users\Pierpaolo\Documents\MATLAB\Microstati\Dati\Results\SubjectFit\Mappe\PazientiRight\F1\4templates';
singleMapRight ='C:\Users\Pierpaolo\Documents\MATLAB\Microstati\Dati\Results\SubjectFit\Mappe\PazientiRight\F2\4templates';

%controlsMapsF1=dir(pathControlsMapsF1);
controlsMaps=dir(singleMapControl );
controlsMaps(1)=[];
controlsMaps(1)=[];
%leftMapsF1=dir(pathLeftMapsF1);
leftMaps=dir(singleMapLeft);
leftMaps(1)=[];
leftMaps(1)=[];
%rightMapsF1=dir(pathRightMapsF1);
rightMaps=dir(singleMapRight);
rightMaps(1)=[];
rightMaps(1)=[];
EEG.chanlocs = readlocs('19ch2.xyz');

%carico i template per gruppo
templatesControlli=load('C:\Users\Pierpaolo\Documents\MATLAB\Microstati\Dati\MappeMedie\groupTemplatesControl.ep');
templatesLeft=load('C:\Users\Pierpaolo\Documents\MATLAB\Microstati\Dati\MappeMedie\groupTemplates4LeftOrd.ep');
templatesRight=load('C:\Users\Pierpaolo\Documents\MATLAB\Microstati\Dati\MappeMedie\groupTemplatesRight.ep');


%         subplot(2,2,4);
%         topoplot(templatesLeft(4,:), EEG.chanlocs)
%ordino le mappe per global dissimilarity
%for i= 3:1%size(leftMaps,1)
clear y
i=30
y = load(strcat('C:\Users\Pierpaolo\Documents\MATLAB\Microstati\Dati\Results\SubjectFit\Mappe\PazientiLeft\F2\4templates\',leftMaps(i).name));
figure
        subplot(2,2,1);topoplot(templatesLeft(1,:), EEG.chanlocs);title('A');
        subplot(2,2,2);topoplot(templatesLeft(2,:), EEG.chanlocs);title('B');
        subplot(2,2,3);topoplot(templatesLeft(3,:), EEG.chanlocs);title('C');    
        subplot(2,2,4);
        topoplot(templatesLeft(4,:), EEG.chanlocs);title('D');
figure
        subplot(2,2,1);topoplot(y(1,:), EEG.chanlocs);title('1');
        subplot(2,2,2);topoplot(y(2,:), EEG.chanlocs);title('2')
        subplot(2,2,3);topoplot(y(3,:), EEG.chanlocs);title('3')
        subplot(2,2,4);topoplot(y(4,:), EEG.chanlocs);title('4')

clear controlliOrdinati
A=y(2,:);
B=y(4,:);
C=y(1,:);
D=y(3,:);
controlliOrdinati=[A;B;C;D];
saveeph(strcat('C:\Users\Pierpaolo\Documents\MATLAB\Microstati\Dati\MappeMedie\subjectWiseOrdinati\left\',leftMaps(i).name),controlliOrdinati)


% figure
%         subplot(2,2,1);
%         topoplot(controlliOrdinati(1,:), EEG.chanlocs)
% 
%         subplot(2,2,2);
%         topoplot(controlliOrdinati(2,:), EEG.chanlocs)
% 
%         subplot(2,2,3);
%         topoplot(controlliOrdinati(3,:), EEG.chanlocs)
%         
%         subplot(2,2,4);
%         topoplot(controlliOrdinati(4,:), EEG.chanlocs)
% % 
templatesRight(1,:)=templatesRight(1,:).*-1
templatesRight(3,:)=templatesRight(3,:).*-1
templatesRight(4,:)=templatesRight(4,:).*-1
    for j=1:4
        
        for jj=1:4
            [diss(j,jj),sc(j,jj)]=computedissandsc(templatesControlli(j,:),templatesRight(jj,:)) ;
        end
        
        %[C,I]=min(diss);
       
        %controlliOrdinati(j,:)=y(I,:);
        %clear diss ind C I
    end
    
    
    
    %check visivo mappe
    f1 = figure;
    movegui(f1,[200,400]);
    subplot(2,2,1);
    topoplot(templatesControlli(1,:), EEG.chanlocs)
    
    subplot(2,2,2);
    topoplot(templatesControlli(2,:), EEG.chanlocs)
    
    subplot(2,2,3);
    topoplot(templatesControlli(3,:), EEG.chanlocs)
    
    subplot(2,2,4);
    topoplot(templatesControlli(4,:), EEG.chanlocs)
%     
    f2 = figure;
    movegui(f2,[800,400]);
    subplot(2,2,1);
    topoplot(templatesRight(1,:), EEG.chanlocs)
    
    subplot(2,2,2);
    topoplot(templatesRight(2,:), EEG.chanlocs)
    
    subplot(2,2,3);
    topoplot(templatesRight(3,:), EEG.chanlocs)
    
    subplot(2,2,4);
    topoplot(templatesRight(4,:), EEG.chanlocs)
    
    % salvo il file per ogni soggetto
    
    
    pause
end

figure
        subplot(2,2,1);
        topoplot(templatesControlli(1,:), EEG.chanlocs)
        title('A')
        subplot(2,2,2);
        topoplot(templatesControlli(2,:), EEG.chanlocs)
title('B')
        subplot(2,2,3);
        topoplot(templatesControlli(3,:), EEG.chanlocs)
 title('C')       
        subplot(2,2,4);
        topoplot(templatesControlli(4,:), EEG.chanlocs)
        title('D')



