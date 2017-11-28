% Calculate repeated measures Anova
%
% [Result,settings] = lab_RepeatedMeasures(data,settings)
%
% data               = array (measures x ratings)
%
% Rm_Anova - Code from Mathworks Fileexchange
%
% Written by F. Hatz 2014

function Result = lab_repeatedmeasures(data,settings)

Result = [];
if ~exist('settings','var')
    settings = [];
end
if ~exist('data','var') | isempty(data)
    resultmode = questdlg('With or without outcomes?','Outcomes?','With','Without','Without');
    settings.collectfiles = 1;
    if strcmp(resultmode,'With')
        [data,header,result,~,settings] = lab_read_statistics(settings,1,0,1,0,1);
    else
        [data,header,~,~,settings] = lab_read_statistics(settings,-1,0,1,0,1);
        result = ones(size(data,1),1);
    end
    if isempty(data)
        return
    end
end
if ~exist('header','var')
    [Filename,Path] = uigetfile('*.xls','Select File to store results');
else
    Filename = header.file;
    Path = header.path;
end
[~,~,~,Filename] = lab_filename(Filename);

if ~isfield(settings,'clustervars') | settings.clustervars == 1
    Prompt = {'Number of Timepoints per measure','clustervars'};
    Formats.type = 'edit';
    Formats.format = 'integer';
    Formats.limits = [2 9999999];
    Formats.size = 30;
    [settings,Cancelled] = inputsdlg(Prompt,'Number of Timepoints per measure',Formats,settings);
    if isempty(settings) | Cancelled == 1 | settings.clustervars <= 1
        return
    else
        pause(0.2);
        settings.numclusters = floor(size(data,2) / settings.clustervars);
    end
end
if ~isfield(settings,'numclusters')
    settings.numclusters = floor(size(data,2) / settings.clustervars);
end

% Create Output-Folder
warning off; %#ok<WNOFF>
mkdir(fullfile(Path,'rmAnova'));
warning on; %#ok<WNON>
Path = fullfile(Path,'rmAnova');

% Reshape data
data = reshape(data,[size(data,1),settings.clustervars,settings.numclusters]);

% Split in groups
if size(result,2) > 1
    disp('Multiple outcomes selected, take only first')
    result = result(:,1);
end
Igroups = unique(result);
Input = cell(size(data,3),length(Igroups));
for j = 1:size(data,3)
    for i = 1:length(Igroups)
        Input{j,i} = data(result == Igroups(i),:,j);
    end
end

% Create names for measures and variables
if isfield(header,'vars')
    varstmp = reshape(header.vars,[size(data,2),size(data,3)]);
    vars = cell(1,size(varstmp,1));
    for i = 1:size(varstmp,1)
        tmp = strfind(varstmp{i,1},'_');
        if ~isempty(tmp)
            vars{1,i} = regexprep(varstmp{i,1}(tmp(end)+1:end),{'*','/','_'},' ');
        else
            vars{1,i} = regexprep(varstmp{i,1},{'*','/','_'},' ');
        end
    end
    Result.Timepoints = vars;
    clearvars vars
    
    vars = cell(1,size(varstmp,2));
    for i = 1:size(varstmp,2)
        tmp = strfind(varstmp{1,i},'_');
        if ~isempty(tmp)
            vars{1,i} = regexprep(varstmp{1,i}(1:tmp(end)-1),{'*','/','_'},' ');
        else
            vars{1,i} = regexprep(varstmp{1,i},{'*','/','_'},' ');
        end
    end
    Measures = vars;
    clearvars vars varstmp i
else
    for i = 1:size(data,3)
        Result.Measure{1,i} = ['Measure' num2str(i)];
    end
    for i = 1:size(data,2)
        Result.Variables{1,i} = ['Rating' num2str(i)];
    end
    clearvars i
end

NumMeasures = size(Input,1);
for Nm = 1:NumMeasures
    [p,xlsout] = anova_rm(Input(Nm,:),'off');
    Result.(Measures{Nm}).p = p;
    Result.(Measures{Nm}).table = xlsout;
    lab_write_xls(fullfile(Path,[Filename '.xlsx']),xlsout,Measures{Nm});
end

end