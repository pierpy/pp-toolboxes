clear all
close all
clc

path=strcat('C:\Users\Pierpaolo\Documents\MATLAB\MicrostatiBetaMotorio\matant_microstati\tutti');
cd(path);
files=dir('*.txt');


    for j = 1:44
%     mov1 = load(strcat(path,'\Mov',num2str(k),num2str(j),'.txt'))
    disp(strcat(path,'\Rest',num2str(j),'.txt'))
    disp(strcat(path,'\Mov1',num2str(j),'.txt'))
    disp(strcat(path,'\Mov2',num2str(j),'.txt'))
    disp(strcat(path,'\Mov3',num2str(j),'.txt'))
    disp(strcat(path,'\Mov4',num2str(j),'.txt'))
    rest = load(strcat(path,'\Rest',num2str(j),'.txt'));
    mov1 = load(strcat(path,'\Mov1',num2str(j),'.txt'));
    mov2 = load(strcat(path,'\Mov2',num2str(j),'.txt'));
    mov3 = load(strcat(path,'\Mov3',num2str(j),'.txt'));
    mov4 = load(strcat(path,'\Mov4',num2str(j),'.txt'));
    
    tot = [rest;mov1;mov2;mov3;mov4];
    
    fid = fopen(strcat('trial',num2str(j),'.txt'),'wt');
        for ii = 1:size(tot,1)
            fprintf(fid,'%10.4f\t',tot(ii,:));
            fprintf(fid,'\n');
        end
        
    fclose(fid)
    
    clear rest mov1 mov2 mov3 mov4 tot
    end

