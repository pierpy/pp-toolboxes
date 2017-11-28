clc
clear all
close all

leftAG_baseline = load('C:\Users\Pierpaolo\Documents\MATLAB\MicrostatiStimolazione\mean_templates\Left-AG\4mean_baseline.txt');
leftAG_time1 = load('C:\Users\Pierpaolo\Documents\MATLAB\MicrostatiStimolazione\mean_templates\Left-AG\4mean_time_1.txt');
leftAG_time2 = load('C:\Users\Pierpaolo\Documents\MATLAB\MicrostatiStimolazione\mean_templates\Left-AG\4mean_time_2.txt');

leftIPS_baseline = load('C:\Users\Pierpaolo\Documents\MATLAB\MicrostatiStimolazione\mean_templates\Left-IPS\4mean_baseline.txt');
leftIPS_time1 = load('C:\Users\Pierpaolo\Documents\MATLAB\MicrostatiStimolazione\mean_templates\Left-IPS\4mean_time_1.txt');
leftIPS_time2 = load('C:\Users\Pierpaolo\Documents\MATLAB\MicrostatiStimolazione\mean_templates\Left-IPS\4mean_time_2.txt');

rightAG_baseline = load('C:\Users\Pierpaolo\Documents\MATLAB\MicrostatiStimolazione\mean_templates\Right-AG\4mean_baseline.txt');
rightAG_time1 = load('C:\Users\Pierpaolo\Documents\MATLAB\MicrostatiStimolazione\mean_templates\Right-AG\4mean_time_1.txt');
rightAG_time2 = load('C:\Users\Pierpaolo\Documents\MATLAB\MicrostatiStimolazione\mean_templates\Right-AG\4mean_time_2.txt');

rightIPS_baseline = load('C:\Users\Pierpaolo\Documents\MATLAB\MicrostatiStimolazione\mean_templates\Right-IPS\4mean_baseline.txt');
rightIPS_time1 = load('C:\Users\Pierpaolo\Documents\MATLAB\MicrostatiStimolazione\mean_templates\Right-IPS\4mean_time_1.txt');
rightIPS_time2 = load('C:\Users\Pierpaolo\Documents\MATLAB\MicrostatiStimolazione\mean_templates\Right-IPS\4mean_time_2.txt');

sham_baseline = load('C:\Users\Pierpaolo\Documents\MATLAB\MicrostatiStimolazione\mean_templates\Sham\4mean_baseline.txt');
sham_time1 = load('C:\Users\Pierpaolo\Documents\MATLAB\MicrostatiStimolazione\mean_templates\Sham\4mean_time_1.txt');
sham_time2 = load('C:\Users\Pierpaolo\Documents\MATLAB\MicrostatiStimolazione\mean_templates\Sham\4mean_time_2.txt');

EEG.chanlocs = readlocs('nuovo_lay.xyz');

figure
subplot(1,4,1);
topoplot(leftAG_baseline(1,:), EEG.chanlocs)
title('\bf1');
subplot(1,4,2);
topoplot(leftAG_baseline(2,:), EEG.chanlocs)
title('\bf2')
subplot(1,4,3);
topoplot(leftAG_baseline(3,:), EEG.chanlocs)
title('\bf3')
subplot(1,4,4);
topoplot(leftAG_baseline(4,:), EEG.chanlocs)
title('\bf4')

ha1 = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
text(0.5, 1,'\bf leftAG baseline','HorizontalAlignment','center','VerticalAlignment', 'top')

figure
subplot(1,4,1);
topoplot(leftAG_time1(1,:), EEG.chanlocs)
title('\bf1');
subplot(1,4,2);
topoplot(leftAG_time1(2,:), EEG.chanlocs)
title('\bf2')
subplot(1,4,3);
topoplot(leftAG_time1(3,:), EEG.chanlocs)
title('\bf3')
subplot(1,4,4);
topoplot(leftAG_time1(4,:), EEG.chanlocs)
title('\bf4')
ha1 = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
text(0.5, 1,'\bf leftAG time1','HorizontalAlignment','center','VerticalAlignment', 'top')


figure
subplot(1,4,1);
topoplot(leftAG_time2(1,:), EEG.chanlocs)
title('\bf1');
subplot(1,4,2);
topoplot(leftAG_time2(2,:), EEG.chanlocs)
title('\bf2');
subplot(1,4,3);
topoplot(leftAG_time2(3,:), EEG.chanlocs)
title('\bf3');
subplot(1,4,4);
topoplot(leftAG_time2(4,:), EEG.chanlocs)
title('\bf4');
ha4 = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
text(0.5, 1,'\bf leftAG time2','HorizontalAlignment','center','VerticalAlignment', 'top')


figure
subplot(1,4,1);
topoplot(leftIPS_baseline(1,:), EEG.chanlocs)
title('\bf1');
subplot(1,4,2);
topoplot(leftIPS_baseline(2,:), EEG.chanlocs)
title('\bf2');
subplot(1,4,3);
topoplot(leftIPS_baseline(3,:), EEG.chanlocs)
title('\bf3');
subplot(1,4,4);
topoplot(leftIPS_baseline(4,:), EEG.chanlocs)
title('\bf4');
ha4 = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
text(0.5, 1,'\bf leftIPS baseline','HorizontalAlignment','center','VerticalAlignment', 'top')

figure
subplot(1,4,1);
topoplot(leftIPS_time1(1,:), EEG.chanlocs)
title('\bf1');
subplot(1,4,2);
topoplot(leftIPS_time1(2,:), EEG.chanlocs)
title('\bf2');
subplot(1,4,3);
topoplot(leftIPS_time1(3,:), EEG.chanlocs)
title('\bf3');
subplot(1,4,4);
topoplot(leftIPS_time1(4,:), EEG.chanlocs)
title('\bf4');
ha4 = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
text(0.5, 1,'\bf leftIPS time1','HorizontalAlignment','center','VerticalAlignment', 'top')

figure
subplot(1,4,1);
topoplot(leftIPS_time2(1,:), EEG.chanlocs)
title('\bf1');
subplot(1,4,2);
topoplot(leftIPS_time2(2,:), EEG.chanlocs)
title('\bf2');
subplot(1,4,3);
topoplot(leftIPS_time2(3,:), EEG.chanlocs)
title('\bf3');
subplot(1,4,4);
topoplot(leftIPS_time2(4,:), EEG.chanlocs)
title('\bf4');
ha4 = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
text(0.5, 1,'\bf leftIPS time2','HorizontalAlignment','center','VerticalAlignment', 'top')

figure
subplot(1,4,1);
topoplot(rightAG_baseline(1,:), EEG.chanlocs)
title('\bf1');
subplot(1,4,2);
topoplot(rightAG_baseline(2,:), EEG.chanlocs)
title('\bf2');
subplot(1,4,3);
topoplot(rightAG_baseline(3,:), EEG.chanlocs)
title('\bf3');
subplot(1,4,4);
topoplot(rightAG_baseline(4,:), EEG.chanlocs)
title('\bf4');
ha4 = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
text(0.5, 1,'\bf rightAG baseline','HorizontalAlignment','center','VerticalAlignment', 'top')

figure
subplot(1,4,1);
topoplot(rightAG_time1(1,:), EEG.chanlocs)
title('\bf1');
subplot(1,4,2);
topoplot(rightAG_time1(2,:), EEG.chanlocs)
title('\bf2');
subplot(1,4,3);
topoplot(rightAG_time1(3,:), EEG.chanlocs)
title('\bf3');
subplot(1,4,4);
topoplot(rightAG_time1(4,:), EEG.chanlocs)
title('\bf4');
ha4 = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
text(0.5, 1,'\bf rightAG time1','HorizontalAlignment','center','VerticalAlignment', 'top')

figure
subplot(1,4,1);
topoplot(rightAG_time2(1,:), EEG.chanlocs)
title('\bf1');
subplot(1,4,2);
topoplot(rightAG_time2(2,:), EEG.chanlocs)
title('\bf2');
subplot(1,4,3);
topoplot(rightAG_time2(3,:), EEG.chanlocs)
title('\bf3');
subplot(1,4,4);
topoplot(rightAG_time2(4,:), EEG.chanlocs)
title('\bf4');
ha4 = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
text(0.5, 1,'\bf rightAG time2','HorizontalAlignment','center','VerticalAlignment', 'top')

figure
subplot(1,4,1);
topoplot(rightIPS_baseline(1,:), EEG.chanlocs)
title('\bf1');
subplot(1,4,2);
topoplot(rightIPS_baseline(2,:), EEG.chanlocs)
title('\bf2');
subplot(1,4,3);
topoplot(rightIPS_baseline(3,:), EEG.chanlocs)
title('\bf3');
subplot(1,4,4);
topoplot(rightIPS_baseline(4,:), EEG.chanlocs)
title('\bf4');
ha4 = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
text(0.5, 1,'\bf rightIPS baseline','HorizontalAlignment','center','VerticalAlignment', 'top')

figure
subplot(1,4,1);
topoplot(rightIPS_time1(1,:), EEG.chanlocs)
title('\bf1');
subplot(1,4,2);
topoplot(rightIPS_time1(2,:), EEG.chanlocs)
title('\bf2');
subplot(1,4,3);
topoplot(rightIPS_time1(3,:), EEG.chanlocs)
title('\bf3');
subplot(1,4,4);
topoplot(rightIPS_time1(4,:), EEG.chanlocs)
title('\bf4');
ha4 = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
text(0.5, 1,'\bf rightIPS time1','HorizontalAlignment','center','VerticalAlignment', 'top')

figure
subplot(1,4,1);
topoplot(rightIPS_time2(1,:), EEG.chanlocs)
title('\bf1');
subplot(1,4,2);
topoplot(rightIPS_time2(2,:), EEG.chanlocs)
title('\bf2');
subplot(1,4,3);
topoplot(rightIPS_time2(3,:), EEG.chanlocs)
title('\bf3');
subplot(1,4,4);
topoplot(rightIPS_time2(4,:), EEG.chanlocs)
title('\bf4');
ha4 = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
text(0.5, 1,'\bf rightIPS time2','HorizontalAlignment','center','VerticalAlignment', 'top')

figure
subplot(1,4,1);
topoplot(sham_baseline(1,:), EEG.chanlocs)
title('\bf1');
subplot(1,4,2);
topoplot(sham_baseline(2,:), EEG.chanlocs)
title('\bf2');
subplot(1,4,3);
topoplot(sham_baseline(3,:), EEG.chanlocs)
title('\bf3');
subplot(1,4,4);
topoplot(sham_baseline(4,:), EEG.chanlocs)
title('\bf4');
ha4 = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
text(0.5, 1,'\bf sham baseline','HorizontalAlignment','center','VerticalAlignment', 'top')

figure
subplot(1,4,1);
topoplot(sham_time1(1,:), EEG.chanlocs)
title('\bf1');
subplot(1,4,2);
topoplot(sham_time1(2,:), EEG.chanlocs)
title('\bf2');
subplot(1,4,3);
topoplot(sham_time1(3,:), EEG.chanlocs)
title('\bf3');
subplot(1,4,4);
topoplot(sham_time1(4,:), EEG.chanlocs)
title('\bf4');
ha4 = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
text(0.5, 1,'\bf sham time1','HorizontalAlignment','center','VerticalAlignment', 'top')

figure
subplot(1,4,1);
topoplot(sham_time2(1,:), EEG.chanlocs)
title('\bf1');
subplot(1,4,2);
topoplot(sham_time2(2,:), EEG.chanlocs)
title('\bf2');
subplot(1,4,3);
topoplot(sham_time2(3,:), EEG.chanlocs)
title('\bf3');
subplot(1,4,4);
topoplot(sham_time2(4,:), EEG.chanlocs)
title('\bf4');
ha4 = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
text(0.5, 1,'\bf sham time2','HorizontalAlignment','center','VerticalAlignment', 'top')