% Script for statistics of eeg and meg results
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
% In case of multiple results, analysis is done for every
% result with separated outputs
%
% Optionally data can by transformed before processing, available options
% are: z-value, log, logit, lateralization index (L-R/L+R)
%
% Optionally a selection of outcomes can be used as factors to correct the
% input data for by generalized linear regression model (a generalized
% linear regression model is calculated with selected factors as depending
% variables and input data as outcome. Thereafter, only residuals of input
% data are used for analysis
%
% Statistics: T-test, Anova, Spearman, Kendall, Pearson, Mann-Whitney-U,
%             Wilcoxon and KruskalWallis
%
% Written by F. Hatz 2012

function Rstat = lab_statistics

Rstat = [];

% Read xls-file
cfg.subfolder = 'BasicStats';
[data,header,result,factors,cfg] = lab_read_statistics(cfg,1,0,0,1,1);
if isempty(data)
    return
end

% Open Logfile
diary(fullfile(header.path,[header.file '_statistics.log']));

% Ask for statistics
disp ('Ask for statistics method')
for i = 1:size(result,2)
    settings(i).Outcome = regexprep(header.result{1,i},'_',' '); %#ok<AGROW>
    Formats{1} = 'text';
    if length(unique(result(:,i))) == 2
        settings(i).Statistics = 'T-test'; %#ok<AGROW>
        settings(i).Write_Mean_Std = true; %#ok<AGROW>
    elseif length(unique(result(:,i))) < 5
        settings(i).Statistics = 'Anova'; %#ok<AGROW>
        settings(i).Write_Mean_Std = true; %#ok<AGROW>
    else
        settings(i).Statistics = 'Spearman'; %#ok<AGROW>
        settings(i).Write_Mean_Std = false; %#ok<AGROW>
    end
    settings(i).Test_Normal_Distribution = false; %#ok<AGROW>
    Formats{2} =  {'T-test','Paired T-test','Anova','Mann-Whitney-U','Wilcoxon','KruskalWallis','Spearman','Kendall','Pearson'};
end
[settings,Cancelled] = inputsdlg(settings,'Statistics',Formats);
if Cancelled == 1
    return
else
    doNormalDistr = false;
    for i = 1:length(settings)
        cfg.method{i,1} = settings(i).Statistics;
        cfg.permutations(i) = 1;
        cfg.write_meanstd(i) = settings(i).Write_Mean_Std;
        if settings(i).Test_Normal_Distribution == true
            doNormalDistr = true;
        end
    end
end
pause(0.2);

% Do statistics (using permutation script with 1 permutation)
[Rstat,skipprocessing] = lab_permutation_calc(data,result,cfg,factors);
if skipprocessing == 1
    return
end
if iscell(Rstat)
    cfg.resultvars = size(Rstat,2);
else
    cfg.resultvars = 1;
    tmp = Rstat;
    clearvars Rstat
    Rstat{1} = tmp;
    clearvars tmp;
end

Rstat = lab_postprocess_data_statistics(Rstat,header,cfg);

% calculate FDR
for Nresult = 1:cfg.resultvars
    [~,~,Rstat{Nresult}.p_fdr] = fdr_bh(Rstat{Nresult}.p);
end

% Do normal distribution
if doNormalDistr == true
    lab_calculate_NormDistribution(data,0.05,header);
end

% Write results of statistics
lab_write_statistics(Rstat,header,cfg);

disp('Finished')
diary off

return