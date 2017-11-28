clear all
close all
clc

pathGroup1='C:\Users\Pierpaolo\Documents\MATLAB\MicrostatiStimolazione\templates\Right-IPS\Time-2\4templates';

group1=dir(pathGroup1);group1(1)=[];group1(1)=[];
EEG.chanlocs = readlocs('Layout_TMS.xyz');

for  k=1:size(group1,1)
      G1=load(strcat(pathGroup1,'\',group1(k).name));
        
        figure;
        for kk=1:size(G1,1)
        subplot(1,4,kk);
        topoplot(G1(kk,:), EEG.chanlocs)
        title(num2str(kk))
        end
        
        button1 = questdlg('invertire mappa 1?');
        if strcmp(button1, 'Yes')
            disp('inverto mappa1');
            G1(1,:)=G1(1,:).*-1;
        end
        button2 = questdlg('invertire mappa 2?');
        if strcmp(button2, 'Yes')
            disp('inverto mappa2');
             G1(2,:)=G1(2,:).*-1;
        end
        button3 = questdlg('invertire mappa 3?');
        if strcmp(button3, 'Yes')
            disp('inverto mappa3');
            G1(3,:)=G1(3,:).*-1;
        end
        button4 = questdlg('invertire mappa 4?');
        if strcmp(button4, 'Yes')
            disp('inverto mappa4');
            G1(4,:)=G1(4,:).*-1;
        end
         
        close all
        figure;
        for kk=1:size(G1,1)
        subplot(1,4,kk);
        topoplot(G1(kk,:), EEG.chanlocs)
        title(num2str(kk))
        end
        
        answer = inputdlg('inserisci sequenza');
        sequenza = str2num(answer{1});
        a = sequenza(1);
        b = sequenza(2);
        c = sequenza(3);
        d = sequenza(4);
        
        if exist('A_sham_t1', 'var') && exist('B_sham_t1', 'var') && exist('C_sham_t1', 'var') && exist('D_sham_t1', 'var')
            A_sham_t1 = vertcat(A_sham_t1, G1(a,:));
            B_sham_t1 = vertcat(B_sham_t1, G1(b,:));
            C_sham_t1 = vertcat(C_sham_t1, G1(c,:));
            D_sham_t1 = vertcat(D_sham_t1, G1(d,:));
        else
            A_sham_t1 = G1(a,:);
            B_sham_t1 = G1(b,:);
            C_sham_t1 = G1(c,:);
            D_sham_t1 = G1(d,:);
        end
        
%         figure;
%         for kk=1:size(A_sham_t1aseline,1)
%             subplot(1,size(A_sham_t1aseline,1),kk);
%             topoplot(A_sham_t1aseline(kk,:), EEG.chanlocs)
%             title(strcat('A_sham_t1aseline',num2str(kk)))
%         end
%         figure;
%         for kk=1:size(B_sham_t1aseline,1)
%             subplot(1,size(B_sham_t1aseline,1),kk);
%             topoplot(B_sham_t1aseline(kk,:), EEG.chanlocs)
%             title(strcat('B_sham_t1aseline',num2str(kk)))
%         end
%         figure;
%         for kk=1:size(C_sham_t1aseline,1)
%             subplot(1,size(C_sham_t1aseline,1),kk);
%             topoplot(C_sham_t1aseline(kk,:), EEG.chanlocs)
%             title(strcat('C_sham_t1aseline',num2str(kk)))
%         end
%         figure;
%         for kk=1:size(D_sham_t1aseline,1)
%             subplot(1,size(D_sham_t1aseline,1),kk);
%             topoplot(D_sham_t1aseline(kk,:), EEG.chanlocs)
%             title(strcat('D_sham_t1aseline',num2str(kk)))
%         end
        
%         button5 = questdlg('OK ?');
%         if strcmp(button5, 'Yes')
%             clear a b c d G1
%             close all
%         end
        
end
allMaps=zeros(4,27);
%%%solo plottaggio%%%
for  k=1:size(group1,1)
      G1=load(strcat(pathGroup1,'\',group1(k).name));
        allMaps=[allMaps;G1];
%         figure;
%         for kk=1:size(G1,1)
%         subplot(2,8,kk);
%         topoplot(G1(kk,:), EEG.chanlocs)
%         title(num2str(kk))
%         end
        
        
end