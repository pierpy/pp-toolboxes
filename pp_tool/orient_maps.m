function [orientend_maps]=orient_maps(SegmentationResult, nclusters, EEG)
orientend_maps = [];

for  k=1:size(SegmentationResult,2)
      G1=SegmentationResult{1, k}.template{1, nclusters};
        
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
                 
      orientend_maps = [orientend_maps;G1];
      
        clearvars G1
end
end