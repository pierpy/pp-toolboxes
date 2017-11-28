% Function combines calculation of ICC, CronbachAlpha and CoefficentOfVariation
%
% [Result] = lab_test_reliability(data,settings)
%
% written by F. Hatz 2014

function [Result] = lab_test_reliability(data,settings)

if ~exist('settings','var')
    settings = [];
end
if ~exist('data','var') | isempty(data)
    [data,header,~,~,settings] = lab_read_statistics([],-1,0,1,1,1);
    if isempty(data)
        Result = [];
        return
    end
end

if ~isfield(settings,'clustervars') | settings.clustervars == 1
    Prompt = {'Number of Judgements per measure','clustervars'};
    Formats.type = 'edit';
    Formats.format = 'integer';
    Formats.limits = [2 9999999];
    Formats.size = 30;
    [settings,Cancelled] = inputsdlg(Prompt,'Number of Judgements per measure',Formats,settings);
    if isempty(settings) | Cancelled == 1 | settings.clustervars <= 1
        Result = [];
        return
    else
        pause(0.2);
    end
end
if ~isfield(settings,'numclusters')
    settings.numclusters = floor(size(data,2) / settings.clustervars);
end

% Create Output-Folder
warning off; %#ok<WNOFF>
mkdir(fullfile(header.path,'TestReliability'));
warning on; %#ok<WNON>
settings.path = fullfile(header.path,'TestReliability');

% Reshape data and correct for NaN
data = reshape(data,[size(data,1),settings.clustervars,settings.numclusters]);
exclude = [];
for i = 1:size(data,3)
    exclude = union(exclude,find(sum(isnan(data(:,:,i)),2)>0));
end
if size(exclude,1) > 1
    exclude = exclude';
end
if ~isempty(exclude)
    if ~isfield(settings,'button')
        settings.button = questdlg('NaN-values: Exclude or Replace by zero','NaN-values','Exclude','Replace','Replace');
    end
    if strcmp(settings.button,'Exclude')
        disp(['Exclude Subjects/Trials with NaN-values (' num2str(exclude) ')'])
        include = 1:size(data,1);
        include = setdiff(include,exclude);
        data = data(include,:,:);
    else
        disp('Replace NaN-values by 0')
        data(isnan(data)) = 0;
    end
end

Ntrials = size(data,1);
for Nvar = 1:size(data,3)
    datatmp = data(:,:,Nvar);
    for i = 2:Ntrials
        Result.CronbachAlpha(Nvar,i) = CronbachAlpha(datatmp(1:i,:));
    end
    fig1 = figure('Color',[1 1 1],'Visible','off');
    plot(Result.CronbachAlpha(Nvar,:));
    title('CronbachAlpha');
    lab_print_figure(fullfile(settings.path,[header.file '_' header.measures{Nvar} '_CronbachAlpha.tif']),fig1);
    close(fig1);
    
    for i = 2:Ntrials
        Result.ICC(Nvar,i) = ICC(datatmp(1:i,:),'A-1',0.05);
    end
    fig2 = figure('Color',[1 1 1],'Visible','off');
    plot(Result.ICC(Nvar,:));
    title('ICC');
    lab_print_figure(fullfile(settings.path,[header.file '_' header.measures{Nvar} '_ICC.tif']),fig2);
    close(fig2);
    
    for i = 2:Ntrials
        Result.CoV(Nvar,i) = std(datatmp(:)) ./ mean(datatmp(:));
    end
    fig3 = figure('Color',[1 1 1],'Visible','off');
    plot(Result.CoV(Nvar,:));
    title('Coefficient of Variance');
    lab_print_figure(fullfile(settings.path,[header.file '_' header.measures{Nvar} '_CoV.tif']),fig2);
    close(fig3);
end

lab_write_xls(fullfile(settings.path,[header.file '_CronbachAlpha.xls']),num2cell(Result.CronbachAlpha));
lab_write_xls(fullfile(settings.path,[header.file '_ICC.xls']),num2cell(Result.ICC));
lab_write_xls(fullfile(settings.path,[header.file '_CoV.xls']),num2cell(Result.CoV));