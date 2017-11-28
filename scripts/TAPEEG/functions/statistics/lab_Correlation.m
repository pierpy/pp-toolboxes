% Calculate Pearson-correlation
%
% [Result,settings] = lab_Correlation(data,settings)
%
% data = array (subjects x measures)
% settings (optional)
%
% Written by F. Hatz 2014

function [Result,settings] = lab_Correlation(data,settings)

if ~exist('settings','var')
    settings = [];
end
if ~exist('data','var') | isempty(data)
    settings.selectmode = 'Multiple';
    settings.collectfiles = 1;
    [data,header,~,~,settings] = lab_read_statistics(settings,-1,0,1,1,1);
    if isempty(data)
        Result = [];
        return
    end
end
if ~isfield(settings,'clustervars') | settings.clustervars <= 1
    disp('Abort, select at least two files to compare')
    return
end
if mod(size(data,2),settings.clustervars) ~= 0
    disp('Abort, mismatch in data')
    return
end

% get settings
if ~isfield(settings,'type') | isempty(settings.type)
    [settings,Cancelled] = lab_set_correlation(settings);
    if Cancelled == 1
        return
    end
end

if ~isfield(settings,'type') | (~strcmp(settings.type,'Pearson') & ~strcmp(settings.type,'Kendall') & ~strcmp(settings.type,'Spearman'))
    settings.type = 'Pearson';
end
if isfield(settings,'excludeNaN') & settings.excludeNaN == true
    ExcludeNaN = 'complete';
else
    ExcludeNaN = 'all';
end

% prepare data for calculation
data = reshape(data,[size(data,1) settings.clustervars size(data,2)/settings.clustervars]);
data = permute(data,[1 3 4 2]);
if isfield(settings,'clustervars2') & settings.clustervars2 > 1 & mod(size(data,2),settings.clustervars2) == 0
    data = reshape(data,[size(data,1) settings.clustervars2 size(data,2)/settings.clustervars2 size(data,4)]);
    data = permute(data,[1 3 2 4]);
end

% Calculate correlation for meaures (subjects are correlated)
Ndata = size(data,4);
Ntrials = size(data,3);
Result = [];
for T = 1:Ndata-1
    for T2 = T+1:Ndata
        Result(end+1).name = [header.variables{T} ' vs. ' header.variables{T2}];
        for N = 1:Ntrials
            [tmp1,tmp2] = corr(data(:,:,N,T),data(:,:,N,T2),'type',settings.type,'rows',ExcludeNaN);
            Result(end).Corr(:,N) = diag(tmp1);
            Result(end).Pval(:,N) = diag(tmp2);
        end
    end
end

% Calculate correlation for subjects (measures and variables are correlated)
dataS = reshape(data,[size(data,1) size(data,2)*size(data,3) size(data,4)]);
dataS = permute(dataS,[2 1 3]);
Ndata = size(dataS,3);
ResultS = [];
for T = 1:Ndata-1
    for T2 = T+1:Ndata
        ResultS(end+1).name = [header.variables{T} ' vs. ' header.variables{T2}];
        [tmp1,tmp2] = corr(dataS(:,:,T),dataS(:,:,T2),'type',settings.type,'rows',ExcludeNaN);
        ResultS(end).Corr = diag(tmp1);
        ResultS(end).Pval = diag(tmp2);
    end
end

% Create Output-Folder
warning off; %#ok<WNOFF>
mkdir(fullfile(header.path,'Corr'));
warning on; %#ok<WNON>
settings.path = fullfile(header.path,'Corr');

% write results
[clustervars,vars] = getstructure(header.measures);
if clustervars > 1 & mod(length(header.measures),clustervars) == 0
    xlsout = {};
    for i = 1:length(Result)
        tmpout = cat(1,{Result(i).name;''},vars);
        tmpoutM = {'Mean'};
        if isfield(header,'variables2')
            tmpoutM(2,1:length(header.variables2)) = header.variables2;
        else
            tmpoutM(2,1) = {''};
        end
        for j = 1:length(vars)
            tmpoutM(2+j,:) = num2cell(mean(Result(i).Corr((j-1)*clustervars+1:j*clustervars,:),1));
        end
        tmpout = cat(2,tmpout,tmpoutM);
        tmpout(1,end+1) = {''};
        tmpoutS = {'Std'};
        if isfield(header,'variables2')
            tmpoutS(2,1:length(header.variables2)) = header.variables2;
        else
            tmpoutS(2,1) = {''};
        end
        for j = 1:length(vars)
            tmpoutS(2+j,:) = num2cell(std(Result(i).Corr((j-1)*clustervars+1:j*clustervars,:),[],1));
        end
        tmpout = cat(2,tmpout,tmpoutS);
        tmpout(1,end+1) = {''};
        tmpoutP = {'P-value'};
        if isfield(header,'variables2')
            tmpoutP(2,1:length(header.variables2)) = header.variables2;
        else
            tmpoutP(2,1) = {''};
        end
        for j = 1:length(vars)
            tmpoutP(2+j,:) = num2cell(mean(Result(i).Pval((j-1)*clustervars+1:j*clustervars,:),1));
        end
        tmpout = cat(2,tmpout,tmpoutP);
        xlsout = cat(1,xlsout,tmpout);
        if i ~= length(Result)
            xlsout(end+1,1) = {''}; %#ok<AGROW>
        end
    end
    if isfield(header,'file')
        lab_write_xls(fullfile(settings.path,[header.file '_CorrCluster.xlsx']),xlsout);
    end
end

% write measures correlations
xlsout = cat(1,{'';''},header.measures);
for i = 1:length(Result)
    tmpout = {Result(i).name};
    if isfield(header,'variables2')
        tmpout(2,1:length(header.variables2)) = header.variables2;
    else
        tmpout(2,1) = {''};
    end
    tmpout = cat(1,tmpout,num2cell(Result(i).Corr));
    xlsout = cat(2,xlsout,tmpout);
end
if isfield(header,'file')
    lab_write_xls(fullfile(settings.path,[header.file '_Corr.xlsx']),xlsout);
end

% write measures p-values
xlsout = cat(1,{'';''},header.measures);
for i = 1:length(Result)
    tmpout = {Result(i).name};
    if isfield(header,'variables2')
        tmpout(2,1:length(header.variables2)) = header.variables2;
    else
        tmpout(2,1) = {''};
    end
    tmpout = cat(1,tmpout,num2cell(Result(i).Pval));
    xlsout = cat(2,xlsout,tmpout);
end
if isfield(header,'file')
    lab_write_xls(fullfile(settings.path,[header.file '_Pvalue.xlsx']),xlsout);
end

% write subjects correlations
xlsout = cat(1,{'Correlation'},header.subjects);
for i = 1:length(ResultS)
    xlsout{1,i+1} = {ResultS(i).name};
    xlsout(2:end,i+1) = num2cell(ResultS(i).Corr);
end
if isfield(header,'file')
    lab_write_xls(fullfile(settings.path,[header.file '_Subjects_Corr.xlsx']),xlsout);
end

% write subjects p-values
xlsout = cat(1,{'Pvalue'},header.subjects);
for i = 1:length(ResultS)
    xlsout{1,i+1} = {ResultS(i).name};
    xlsout(2:end,i+1) = num2cell(ResultS(i).Pval);
end
if isfield(header,'file')
    lab_write_xls(fullfile(settings.path,[header.file '_Subjects_Pval.xlsx']),xlsout);
end

end

function [clustervars,vars] = getstructure(vars)
    for i = 1:length(vars)
        tmp = strfind(vars{i},'_');
        if ~isempty(tmp)
            tmp2{i,1} = vars{i}(1:tmp(end)-1); %#ok<AGROW>
        else
            tmp2{i,1} = vars{i}; %#ok<AGROW>
        end
    end
    [~,m,~]=unique(tmp2,'last');
    m = sort(m);
    n = mod(m,m(1));
    if length(m) == 1 | (m(1) > 1 & (m(2)-m(1)) == m(1))
        clustervars = m(1);
    else
        clustervars = 1;
    end
    clearvars tmp tmp2 m n i
    if clustervars > 1
        for i = 1:floor(length(vars)/clustervars)
            tmp = strfind(vars{(i-1)*clustervars+1},'_');
            if ~isempty(tmp)
                tmp2{i,1} = vars{(i-1)*clustervars+1}(1:tmp(end)-1);
            else
                tmp2{i,1} = vars{(i-1)*clustervars+1};
            end
        end
        vars = tmp2;
    end
end
