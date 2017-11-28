
clear all
close all
clc

files_path='C:\Users\Pierpaolo\Documents\MATLAB\MicrostatiBetaMotorio\Dati_Microstati_MOTORIO';

prompt = 'cartella nome soggetto:';
dlg_title = 'Input';
str = input(prompt,'s');
files_path1=strcat(files_path,'\',str);
cd(files_path1)

load(strcat(str(1:6),'_emg.mat'));



 
        bch = 1:129;
        bch(gch)=[];
        bch(find(bch==48))=[];
        bch(find(bch==119))=[];
        bch(find(bch==125))=[];
        bch(find(bch==126))=[];
        bch(find(bch==127))=[];
        bch(find(bch==128))=[];

        ind = find((bch)> 48);
        bch(ind)=bch(ind)-1;
        ind1 = find((bch)>118);
        bch(ind1)=bch(ind1)-1;
        ind2=find((bch)>122);
        bch(ind2)= bch(ind2)-4;
bch
