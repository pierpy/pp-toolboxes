% Script for permutation of eeg and meg results
% Input data: xls-files or output of Brainwave (Stam, Amsterdam)
%             You can use separated files for different groups of
%             patients / subjects, define groups manually or use a
%             corresponding xls-file with outcomes (subject names
%             must be equal)
% xls-file format: 1.row:     header with subjects
%                  1. column: names of variables
%                              - order: var1ch1,var1ch2,var1ch3...var2ch1..
%                                (var=variable ch=channel)
% Optionally the last rows can contain outcomes (recognized in case of
% multichannel data)
%
% In case of multiple results, permutation analysis is done for every
% result with separated outputs
%
% Optionally data can by transformed before processing, available options
% are: z-value, log, logit, lateralization index (L-R/L+R)
%
% Optionally a selection of outcomes can be used as factors to correct the
% input data for by generalized linear regression model (a generalized
% linear regression model is calculated with selected factors as depending
% variables and input data as outcome. Thereafter, only residuals of input
% data are used for permutation analysis
%
% Permutation Statistics: T-test, Anova, Spearman, Kendall, Pearson, 
%                         Mann-Whitney-U, Wilcoxon and KruskalWallis
%
% After permutation results can be plotted in signal (eeglab routines) or
% source space
% For source space plotting the file 'MNIbrain.mat' is used (see
% documentation)
%
% Written by F. Hatz 2012

function Rstat = lab_permutation

Rstat = [];
    
% Read xls-file
cfg.subfolder = 'Permutation';
[data,header,result,factors,cfg] = lab_read_statistics(cfg,1,1,1,1,1);
if isempty(data)
    return
end

% Open Logfile
diary(fullfile(header.path,[header.file '_permutation.log']));

% Ask for statistics
disp ('Ask for statistics method')

% if isempty(cfg.clustervars) | cfg.clustervars > 1
%     Permutations = 10000;
% else
%     Permutations = 1;
% end
Permutations = 10000;

for i = 1:size(result,2)
    settings(i).Measure = regexprep(header.result{1,i},'_',' '); %#ok<AGROW>
    Formats{1} = 'text';
    if length(unique(result(:,i))) == 2
        settings(i).Statistics = 'T-test'; %#ok<AGROW>
        settings(i).Permutations = Permutations; %#ok<AGROW>
        settings(i).Write_Mean_Std = true; %#ok<AGROW>
    elseif length(unique(result(:,i))) < 5
        settings(i).Statistics = 'Anova'; %#ok<AGROW>
        settings(i).Permutations = Permutations; %#ok<AGROW>
        settings(i).Write_Mean_Std = true; %#ok<AGROW>
    else
        settings(i).Statistics = 'Spearman'; %#ok<AGROW>
        settings(i).Permutations = Permutations; %#ok<AGROW>
        settings(i).Write_Mean_Std = false; %#ok<AGROW>
    end
    Formats{2} = {'T-test','Paired T-test','Anova','RM-Anova','Mann-Whitney-U','Wilcoxon','KruskalWallis','Spearman','Kendall','Pearson'};
end
[settings,Cancelled] = inputsdlg(settings,'Permutation',Formats);
if Cancelled == 1
    return
else
    for i = 1:length(settings)
        cfg.method{i,1} = settings(i).Statistics;
        cfg.permutations(i) = settings(i).Permutations;
        cfg.write_meanstd(i) = settings(i).Write_Mean_Std;
    end
    clearvars Permutations settings Cancelled
end
pause(0.2);

% Do permutation
[Rstat,skipprocessing] = lab_permutation_calc(data,result,cfg,factors);
if skipprocessing == 1
    return
end
if iscell(Rstat)
    cfg.resultvars = length(Rstat);
else
    cfg.resultvars = 1;
    tmp = Rstat;
    clearvars Rstat
    Rstat{1} = tmp;
    clearvars tmp;
end
if length(cfg.method) == 1 & strcmp(cfg.method{1,1},'RM-Anova')
    if isfield(header,'variables2') & ~isempty(header.variables2)
        header.variables = header.variables2;
        header.variables2 = header.measures;
        header.measures = cfg.AnovaVars;
        cfg.clustervars = length(header.variables);
        cfg.clustervars2 = length(header.variables2);
        header.vars = {};
        for i = 1:length(header.measures)
            for j = 1:length(header.variables2);
                for m = 1:length(header.variables);
                    header.vars = cat(2,header.vars,{[header.measures{i} '_' header.variables2{j} '_' header.variables{m}]});
                end
            end
        end
    else
        header.variables = header.measures;
        header.measures = cfg.AnovaVars;
        cfg.clustervars = length(header.variables);
        cfg.clustervars2 = 1;
        header.vars = {};
        for i = 1:length(header.measures)
            for m = 1:length(header.variables);
                header.vars = cat(2,header.vars,{[header.measures{i} '_' header.variables{m}]});
            end
        end
    end
end

% Postprocess results of permutation
Rstat = lab_postprocess_data_statistics(Rstat,header,cfg);

% Write results of permutation
lab_write_statistics(Rstat,header,cfg);

if cfg.clustervars == 1 % no plotting of 1 value possible
    disp('Finished')
    diary off
    return
end

disp('Finished')
diary off

return