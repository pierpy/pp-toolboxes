function [ EEG ] = itab_icaclassification( EEG )
%ITAB_ICACLASSIFICATION Summary of this function goes here
    %   Detailed explanation goes here
    disp('IC identification...');

    [nIC, nptIC] = size(EEG.ICA.IC);
    clas = zeros(1, nIC);
%     best_iter = ICA.best_iter.index;
%     clas(ICA.iteration(best_iter).brain_ic) = 1;

    time = EEG.ICA.timeica;
    %-----------------------------------------------
    %------ POWER SPECTRAL DENSITY ESTIMATION  -----
    %-----------------------------------------------

    win = 1024*ceil(4/(EEG.decimationfactor + 1));

    for ix=1:nIC
        [pxx, F] = pwelch(EEG.ICA.IC(ix,:), win,[],[], EEG.srate);
        mspettro(ix,:) = sqrt((EEG.srate/win)*pxx)';
    end
    window = 10;          %expressed in seconds
    window_base = find(time >= time(1) & time < window + time(1));
    flag = 0;
    while flag == 0

        b_ic = [];
        g_ic = [];

        for i = 1: nIC
            sig = EEG.ICA.IC(i,:);
            
            figure;
            subplot(2,5,1:3), plot(time, sig);
            axis tight; title('time course'); xlabel('s');
            subplot(2,2,3), plot(time(window_base), sig(window_base));
            axis tight; title('time course'); xlabel('seconds');
            subplot(2,2,4), plot(F, mspettro(i,:),'r');
            axis([0 min(100, EEG.lfiltfreq) 0 1.1*max(mspettro(i,:))]); title('power spectrum'); xlabel('Hertz');
            subplot(2,5,4:5), topoplot(EEG.ICA.A(:,i), EEG.chanlocs);
%             colormap(Cx_map);
            
            title(strcat('IC no. ',num2str(i)),'FontSize',12);
            %                  print(strcat(dir2, nome_file, '_ic', num2str(i), '_classification.tif'),'-dtiff','-r200');
            
            pause;
            
            if clas(i) > 0
                ButtonName=questdlg('How do you classify this component?', ...
                    'Independent Components classification','Brain signal','Artifact','Brain signal');
                disp(['component no. ' num2str(i) ' - automatic classification : ',ButtonName]);
                
            else
                ButtonName=questdlg('How do you classify this component?', ...
                    'Independent Components classification','Brain signal','Artifact','Artifact');
                disp(['component no. ' num2str(i) ' - automatic classification : ',ButtonName]);
            end
            
            if strcmp(ButtonName,'Brain signal')
                g_ic = [g_ic i];
            else
                b_ic = [b_ic i];
            end
            close all;
        end
        
        clear sig mat;
        
        ButtonName=questdlg('Are you sure that the classification is correct?', ...
            'Independent Components classification','Yes','No','Yes');
        
        if strcmp(ButtonName,'Yes')
            flag=1;
        end
        
    end
    EEG.ICA.gIC = g_ic;
    EEG.ICA.bIC = b_ic;
    
end

