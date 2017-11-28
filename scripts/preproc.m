clc 
clear all
close all

matlab_file_path = 'C:\Users\Pierpaolo\Documents\MATLAB\_Microstati-Prestimulus\Att-AG\angser-att-ag\conc.mat';
chanlocs = 'C:\Users\Pierpaolo\Documents\MATLAB\_Microstati-Prestimulus\Layout_pre_stimulus_ok.xyz';
nbchan = 28;
pnts = 50;
srate = 256; 

eeg = pop_importdata('data', matlab_file_path,...
                     'dataformat', 'matlab',...
                     'chanlocs', chanlocs,...
                     'nbchan', nbchan,...
                     'pnts', pnts,...
                     'srate',srate );
                 
 %% interpolazione
 eeg_interp = pop_interp(eeg, [7 11 15 27], 'spherical');
