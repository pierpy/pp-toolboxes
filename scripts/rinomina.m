clc
clear all
close all

path='C:\Users\Pierpaolo\Documents\MATLAB\Microstati\DATI_RAGU\right';
cd(path)
files=dir(path);
files(1)=[];files(1)=[];
for i=1:size(files,1)
    nome=files(i).name;
   movefile(nome,strcat(nome(5:10),'_right.asc')) 
end