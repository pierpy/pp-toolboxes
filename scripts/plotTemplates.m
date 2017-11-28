clc
clear all
close all

Templates3Control =load('C:\Users\Pierpaolo\Documents\MATLAB\Microstati\Dati\Results\RESULT LAST\templates\3mappe\group3TemplatesControl.ep');
Templates3Left =load('C:\Users\Pierpaolo\Documents\MATLAB\Microstati\Dati\Results\RESULT LAST\templates\3mappe\groupTemplatesLeft.ep');
Templates3Right =load('C:\Users\Pierpaolo\Documents\MATLAB\Microstati\Dati\Results\RESULT LAST\templates\3mappe\group3TemplatesRight.ep');

Templates4Control =load('C:\Users\Pierpaolo\Documents\MATLAB\Microstati\Dati\Results\RESULT LAST\templates\4mappe\groupTemplatesControl.ep');
Templates4Left =load('C:\Users\Pierpaolo\Documents\MATLAB\Microstati\Dati\Results\RESULT LAST\templates\4mappe\groupTemplates4LeftOrd.ep');
Templates4Right =load('C:\Users\Pierpaolo\Documents\MATLAB\Microstati\Dati\Results\RESULT LAST\templates\4mappe\groupTemplatesRight.ep');

EEG.chanlocs = readlocs('19ch2.xyz');

figure
subplot(2,2,1);
topoplot(Templates3Control(1,:), EEG.chanlocs)
title('\bfA');
subplot(2,2,2);
topoplot(Templates3Control(2,:), EEG.chanlocs)
title('\bfB')
subplot(2,2,3);
topoplot(Templates3Control(3,:), EEG.chanlocs)
title('\bfD')
ha1 = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
text(0.5, 1,'\bf Controlli','HorizontalAlignment','center','VerticalAlignment', 'top')

figure
subplot(2,2,1);
topoplot(Templates3Left(1,:), EEG.chanlocs)
title('\bfA');
subplot(2,2,2);
topoplot(Templates3Left(2,:), EEG.chanlocs)
title('\bfB');
subplot(2,2,3);
topoplot(Templates3Left(3,:), EEG.chanlocs)
title('\bfD');
ha2 = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
text(0.5, 1,'\bf Pazienti Left','HorizontalAlignment','center','VerticalAlignment', 'top')

figure
subplot(2,2,1);
topoplot(Templates3Right(1,:), EEG.chanlocs)
title('\bfA');
subplot(2,2,2);
topoplot(Templates3Right(2,:), EEG.chanlocs)
title('\bfB');
subplot(2,2,3);
topoplot(Templates3Right(3,:), EEG.chanlocs)
title('\bfD');
ha3 = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
text(0.5, 1,'\bf Pazienti Right','HorizontalAlignment','center','VerticalAlignment', 'top')

figure
subplot(2,2,1);
topoplot(Templates4Control(1,:), EEG.chanlocs)
title('\bfA');
subplot(2,2,2);
topoplot(Templates4Control(2,:), EEG.chanlocs)
title('\bfB');
subplot(2,2,3);
topoplot(Templates4Control(3,:), EEG.chanlocs)
title('\bfC');
subplot(2,2,4);
topoplot(Templates4Control(4,:), EEG.chanlocs)
title('\bfD');
ha4 = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
text(0.5, 1,'\bf Controlli','HorizontalAlignment','center','VerticalAlignment', 'top')


figure
subplot(2,2,1);
topoplot(Templates4Left(1,:), EEG.chanlocs)
title('\bfA');
subplot(2,2,2);
topoplot(Templates4Left(2,:), EEG.chanlocs)
title('\bfB');
subplot(2,2,3);
topoplot(Templates4Left(3,:), EEG.chanlocs)
title('\bfC');
subplot(2,2,4);
topoplot(Templates4Left(4,:), EEG.chanlocs)
title('\bfD');
ha5 = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
text(0.5, 1,'\bf Pazienti Left','HorizontalAlignment','center','VerticalAlignment', 'top')

figure
subplot(2,2,1);
topoplot(Templates4Right(1,:), EEG.chanlocs)
title('\bfA');
subplot(2,2,2);
topoplot(Templates4Right(2,:), EEG.chanlocs)
title('\bfB');
subplot(2,2,3);
topoplot(Templates4Right(3,:), EEG.chanlocs)
title('\bfC');
subplot(2,2,4);
topoplot(Templates4Right(4,:), EEG.chanlocs)
title('\bfD');
ha6 = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
text(0.5, 1,'\bf Pazienti Right','HorizontalAlignment','center','VerticalAlignment', 'top')