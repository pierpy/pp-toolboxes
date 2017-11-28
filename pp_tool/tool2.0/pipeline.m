%struttura con informazioni circa l'analisi 
clear all
close all
clc
ff = dir('*.mat');
warning off
confstruct.substract_column_at_start = 1; % detrend si o no
confstruct.method_GFPeak = 'GFPL2';       % gfp calculation method
confstruct.use_gfp_peaks = 1;             % utililla i picchi del gfp o tutto il file
confstruct.setall1 = 0;
confstruct.normalize = 0;
confstruct.tf = 250;
confstruct.fs = 125;
confstruct.nch = 62;
confstruct.minclusters = 12;
confstruct.maxclusters = 12;
confstruct.pmode = 0;
confstruct.algorithm = 1;
confstruct.similarity_measure = 'correlation';
confstruct.debug = 0;
confstruct.segment_min_len = 3;
confstruct.downsample = 1;
RES_parameter = {};
RES_clustering = {};
TOTAL_RES_clustering = {};
TOTAL_RES_parameter = {};
load('C:\Users\Pierpaolo\Documents\MATLAB\Intelligenza\risultati\analisi1_templates_medi\mappeMedieOrdinate.mat');
  maps = b;
for f = 1:1%size(ff)
    load(ff(f,1).name);
    for cond = 1:1%size(DATA,1)
        for trial = 1:1%size(DATA(cond,:),2)
            if ~isempty(DATA{cond,trial}) && size(DATA{cond,trial}.dataB,2) > 4 && sum(sum(abs(DATA{cond, trial}.dataB)))>2
                
                disp(strcat('soggetto', num2str(ff(f,1).name)))
                disp(strcat('condizione', num2str(cond)))
                disp(strcat('trial', num2str(trial)))
                
                eeg_or = downsample(DATA{cond, trial}.dataB',4);
              
                [ eeg, gfp_peaks_indices, gfp_curve ] = modmaps_preprocess (eeg_or, confstruct);

%                 [ NMicr ] = CW2( eeg ,gfp_peaks_indices, confstruct);
% 
%                 maps = NMicr.template{1, 1};
% 
                 [ output ] = compute_mstate_parameters( confstruct, eeg, maps );
                 RES_parameter{cond, trial} = output;
%                 RES_clustering{cond, trial} = NMicr;
%                 
                clear output NMicr
            end
        end
    end
  % TOTAL_RES_clustering_baseline{f,1}=RES_clustering;
    TOTAL_RES_parameter_baseline{f,1}=RES_parameter;
end

% for f = 1:size(ff)
%     load(ff(f,1).name);
%     for cond = 1:size(DATA,1)
%         for trial = 1:size(DATA(cond,:),2)
%             if ~isempty(DATA{cond,trial}) && size(DATA{cond,trial}.dataT,2) > 4 && sum(sum(abs(DATA{cond, trial}.dataT)))>2
%                 
%                 disp(strcat('soggetto', num2str(ff(f,1).name)))
%                 disp(strcat('condizione', num2str(cond)))
%                 disp(strcat('trial', num2str(trial)))
%                 
%                 eeg_or = downsample(DATA{cond, trial}.dataT',4);
%                 
%                 [ eeg, gfp_peaks_indices, gfp_curve ] = modmaps_preprocess (eeg_or, confstruct);
% 
%                 [ NMicr ] = CW2( eeg ,gfp_peaks_indices, confstruct);
% 
% %                 maps = NMicr.template{1, 1};
% % 
% %                 [ output ] = compute_mstate_parameters( confstruct, eeg, maps );
% %                 RES_parameter{cond, trial} = output;
%                 RES_clustering{cond, trial} = NMicr;
% %                 
%                 clear output NMicr
%             end
%         end
%     end
%     TOTAL_RES_clustering_conditions{f,1}=RES_clustering;
%     TOTAL_RES_parameter_conditions{f,1}=RES_parameter;
% end

% 
% mean_dur_accross_trial = [];
% mean_single_var = [];
% for sub = 1:size(TOTAL_RES_parameter_baseline,1)
%     RES_parameter = TOTAL_RES_parameter_baseline{sub,1};
%     for c = 1: size(RES_parameter,1)
%        for i = 1: size(RES_parameter(c,:), 2)
%             if ~isempty(RES_parameter{c,i})
%                 mean_dur_accross_trial(:,i) =  nanmean(RES_parameter{c, i}.mean_dur ,2);
%                 mean_freq_accross_trial(:,i) = nanmean(RES_parameter{c, i}.freq ,2);
%                 mean_cov_accross_trial(:,i) = nanmean(RES_parameter{c,i}.cov,2);
%                 mean_single_var = [mean_single_var; RES_parameter{c,i}.var_ms];   
%             end
%        end
%         mean_dur_cond(:,c) = nanmean(mean_dur_accross_trial,2);
%         mean_freq_cond(:,c) = nanmean(mean_freq_accross_trial, 2);
%         mean_cov(:,c) = nanmean(mean_cov_accross_trial, 2);
%         mean_var(c,:) = nanmean(mean_single_var);  
%     end
%     MEDIE{sub}.dur = mean_dur_cond;
%     MEDIE{sub}.freq = mean_freq_cond;
%     MEDIE{sub}.cov = mean_cov;
%     MEDIE{sub}.var = mean_var;
% end
% % 
% % all1=[];
% all2=[];
% all3=[];
% for k = 2:size(MEDIE,2)
%     all1 = [all1, MEDIE{1,k}.dur(:,1)];
%     all2 = [all2, MEDIE{1, k}.dur(:,2)];
%     all3 = [all3, MEDIE{1, k}.dur(:,3)];
% end
% 
% % 
% for s = 1: size(TOTAL_RES_parameter,1)
%     RES_parameter = TOTAL_RES_parameter{s,1};
%     for c = 1: size(RES_parameter,1)
%         for t = 1:size(RES_parameter, 2)
%          if   ~isempty(RES_parameter{c, t})
%             N = isnan(RES_parameter{c, t}.mean_dur);
%             d = find(N == 1);
%             if ~isempty(d)
%                 disp(strcat('condizione: ', num2str(c)))
%                 disp(strcat('tiral: ', num2str(t)))
%                 disp('NAN')
%             else
%                 disp(strcat('condizione: ', num2str(c)))
%                 disp(strcat('tiral: ', num2str(t)))
%                 disp('OK')
%             end
%             pause
%          end
%        end
%     end
%     disp(strcat('soggetto: ', num2str(s)))
%     pause
% end

