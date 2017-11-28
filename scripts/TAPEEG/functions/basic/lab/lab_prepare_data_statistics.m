% Helper function for lab_read_statistics
%
% written by F. Hatz 2012

function [data,header,result,factors,cfg] = lab_prepare_data_statistics(inputdata,cfg,nofactors,noclusters)

if ~exist('nofactors','var') | isempty(nofactors)
    nofactors = 0;
end
if ~exist('noclusters','var') | isempty(noclusters)
    noclusters = 0;
end

if ~exist('inputdata','var') | isempty(inputdata)
    disp('Input file is invalid')
    data = [];
    result = [];
    factors = [];
    header = [];
    return
else
    inputdata = inputdata';
end

tmp = inputdata{1,1};
if ~isempty(tmp)
    tmp = textscan(tmp,'%s');
    tmp = tmp{1,1};
    if strcmp(tmp{1,1}(1),'C')
        cfg.clustervars = str2num(tmp{1,1}(2:end)); %#ok<ST2NM>
    end
    if size(tmp,1) > 1 & strcmp(tmp{2,1}(1),'R')
        cfg.numresults = str2num(tmp{2,1}(2:end)); %#ok<ST2NM>
    end
end
if noclusters == 1
    cfg.clustervars = 1;
end

if exist('cfg','var') & isfield(cfg,'file')
    header.file = cfg.file;
    header.path = cfg.path;
else
    header.file = [];
    header.path = [];
end
try
    data = cell2mat(inputdata(2:end,2:end));
catch %#ok<CTCH>
    disp('Input file has invalid / non numeric data')
    data = [];
    result = [];
    factors = [];
    header = [];
    return
end

% Ask for settings of input data
disp ('Settings for data')
if ~isfield(cfg,'clustervars')
    cfg.clustervars = [];
end
if ~isfield(cfg,'numresults')
    cfg.numresults = [];
end
if ~isfield(cfg,'calcratio')
    cfg.calcratio = false;
    cfg.excludeM = false;
    cfg.excludeV = false;
    cfg.Latindex = false;
    cfg.zvalues = false;
    cfg.Logdata = false;
    cfg.Logitdata = false;
end

Prompt = {};
Formats = {};

if noclusters == 0
    Prompt(end+1,:) = {'Number of variables per measure', 'clustervars'};
    Formats(end+1,1).type = 'edit';
    Formats(end,1).format = 'integer';
    Formats(end,1).limits = [0 9999];
    Formats(end,1).size = 60;
end

Prompt(end+1,:) = {'Number of outcomes', 'numresults'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'integer';
Formats(end,1).limits = [0 9999];
Formats(end,1).size = 60;

Prompt(end+1,:) = {'Exclude measures' 'excludeM'};
Formats(end+1,1).type = 'check';
Formats(end,1).size = [-1 -1];

Prompt(end+1,:) = {'Exclude variables' 'excludeV'};
Formats(end+1,1).type = 'check';
Formats(end,1).size = [-1 -1];

Prompt(end+1,:) = {'Calculate ratio' 'calcratio'};
Formats(end+1,1).type = 'check';
Formats(end,1).size = [-1 -1];

Prompt(end+1,:) = {'Lateralization index' 'Latindex'};
Formats(end+1,1).type = 'check';
Formats(end,1).size = [-1 -1];

Prompt(end+1,:) = {'Z-values' 'zvalues'};
Formats(end+1,1).type = 'check';
Formats(end,1).size = [-1 -1];

Prompt(end+1,:) = {'Log data' 'Logdata'};
Formats(end+1,1).type = 'check';
Formats(end,1).size = [-1 -1];

Prompt(end+1,:) = {'Logit data' 'Logitdata'};
Formats(end+1,1).type = 'check';
Formats(end,1).size = [-1 -1];

[cfg,Cancelled] = inputsdlg(Prompt,'Transform data',Formats,cfg);
if Cancelled == 1
    data = [];
    result = [];
    factors = [];
    header = [];
    return
end

% Abort if input not matching
if ~isfield(cfg,'clustervars') | isempty(cfg.clustervars) | size(data,2) <= cfg.clustervars
    disp('Abort: Number of electrodes and data-file not matching')
    data = [];
    result = [];
    factors = [];
    header = [];
    return
end

% find results in data and fill header
if ~isfield(cfg,'numresults')
    cfg.numresults = 1;
end
cfg.numclusters = floor((size(data,2)-cfg.numresults)/cfg.clustervars);
resultvars = size(data,2) - (cfg.numclusters * cfg.clustervars);
strlist = inputdata(1,end-resultvars+1:end);
selection = listdlg('PromptString','Select outcomes(s):','SelectionMode','multiple', ...
    'ListString',strlist,'ListSize',[300 400]);
if isempty(selection)
    data = [];
    result = [];
    factors = [];
    header = [];
    return
end
result = data(:,selection + end-resultvars);
header.result = inputdata(1,selection + end-resultvars);
for i = 1:size(header.result,2)
    header.result{1,i} = regexprep(header.result{1,i},{'::',':',';','/','\',','},'_');
end
cfg.resultvars = size(selection,2);
strlist = setdiff(strlist,strlist(selection));
if ~isempty(strlist) & nofactors == 0
    strlist = [cellstr('none') strlist(:)'];
    tmp = 1:resultvars;
    tmp = setdiff(tmp,selection);
    selection = listdlg('PromptString','Select factor(s)','SelectionMode','multiple', ...
        'ListString',strlist,'ListSize',[300 400]);
    if ~isempty(selection) & min(selection) > 1
        factors = data(:,tmp(selection-1) + end-resultvars);
        header.factors = inputdata(1,tmp(selection-1) + end-resultvars);
        cfg.factorvars = size(header.factors,2);
    else
        factors = [];
        header.factors = [];
        cfg.factorvars = 0;
    end
else
    factors = [];
    header.factors = [];
    cfg.factorvars = 0;
end
data = data(:,1:end-resultvars);
header.vars = inputdata(1,2:end-resultvars);
header.subjects = inputdata(2:end,1);
clearvars strlist tmp selection resultvars
for i = 1:cfg.numclusters
    tmp = strfind(header.vars{1,(i-1)*cfg.clustervars+1},'_');
    if cfg.clustervars == 1
        measures{i,1} = header.vars{1,(i-1)*cfg.clustervars+1}; %#ok<AGROW>
    elseif ~isempty(tmp)
        measures{i,1} = header.vars{1,(i-1)*cfg.clustervars+1}(1:tmp(1)-1); %#ok<AGROW>
    else
        measures{i,1} = ['Measure' num2str(i)]; %#ok<AGROW>
    end
    for j = 1:size(measures,1)
        measures{j,1} = regexprep(measures{j,1},':','_'); %#ok<AGROW>
        measures{j,1} = regexprep(measures{j,1},'-','_'); %#ok<AGROW>
        measures{j,1} = regexprep(measures{j,1},filesep,'_'); %#ok<AGROW>
        measures{j,1} = regexprep(measures{j,1},';','.'); %#ok<AGROW>
    end
end
for i = 1:cfg.clustervars
    tmp = strfind(header.vars{1,i},'_');
    if ~isempty(tmp)
        variables{i,1} = header.vars{1,i}(tmp(1)+1:end); %#ok<AGROW>
    else
        variables{i,1} = ['Var' num2str(i)]; %#ok<AGROW>
    end
end
clearvars tmp i

% Calculate ratio
if isfield(cfg,'calcratio') & cfg.calcratio == true
    selection = listdlg('PromptString','Ratio first measure:','SelectionMode','single', ...
        'ListString',measures,'ListSize',[300 400]);
    selection2 = listdlg('PromptString','Ratio second measure:','SelectionMode','single', ...
        'ListString',measures,'ListSize',[300 400]);
    if ~isempty(selection) & ~isempty(selection2)
        data = real(log(data(:,((selection-1)*cfg.clustervars + 1):(selection*cfg.clustervars)) ./ ...
            data(:,((selection2-1)*cfg.clustervars + 1):(selection2*cfg.clustervars))));
    end
    measures = {['Ratio_' measures{selection} '_' measures{selection2}]};
    cfg.numclusters = 1;
    varstmp = [];
    for i = 1:cfg.clustervars
        tmp = strfind(header.vars{1,i},'_');
        if ~isempty(tmp)
            varstmp = [varstmp cellstr(['Ratio_' header.vars{1,i}(tmp(1)+1:end)])]; %#ok<AGROW>
        else
            varstmp = [varstmp cellstr(['Ratio_ch' num2str(i)])]; %#ok<AGROW>
        end
    end
    header.vars = varstmp;
    clearvars varstmp
end
clearvars i j selection selection2

% Exclude measures
if isfield(cfg,'excludeM') & cfg.excludeM == true
    selection = listdlg('PromptString','Select measures:','SelectionMode','multiple', ...
        'ListString',measures,'InitialValue',1:size(measures,1),'ListSize',[300 400]);
    if ~isempty(selection)
        datatmp = [];
        varstmp = [];
        for i = selection
            datatmp = [datatmp data(:,(i-1)*cfg.clustervars+1:i*cfg.clustervars)]; %#ok<AGROW>
            varstmp = [varstmp header.vars(:,(i-1)*cfg.clustervars+1:i*cfg.clustervars)]; %#ok<AGROW>
        end
        data = datatmp;
        header.varsall = header.vars;
        header.vars = varstmp;
        clearvars datatmp
        measures = measures(selection);
        cfg.numclusters = length(selection);
    end
    clearvars selection
end

% Exclude variables
if isfield(cfg,'excludeV') & cfg.excludeV == true
    if cfg.clustervars == 78
        strdefault = [7,8,9,10,14,15,16,17,18,19,20,22,23,30,31,46,47,48,49,53,54,55,56,57,58,59,61,62,69,70];
    else
        strdefault = 1:cfg.clustervars;
    end
    selection = listdlg('PromptString','Select variables to include:','SelectionMode','multiple', ...
        'ListString',variables,'InitialValue',strdefault,'ListSize',[300 400]);
    if ~isempty(selection)
        cfg.includevars = selection;
        cfg.excludevars = setdiff((1:cfg.clustervars),cfg.includevars);
        datatmp = [];
        varstmp = [];
        for i = 1:cfg.numclusters
            datatmp = [datatmp data(:,cfg.includevars + (i-1)*cfg.clustervars)]; %#ok<AGROW>
            varstmp = [varstmp header.vars(:,cfg.includevars + (i-1)*cfg.clustervars)]; %#ok<AGROW>
        end
        data = datatmp;
        cfg.clustervars = size(cfg.includevars,2);
        variables = variables(cfg.includevars,1);
        header.varsall = header.vars;
        header.vars = varstmp;
        clearvars datatmp varstmp
    end
    clearvars selection strdefault
end

% Calculate lateralization index
if isfield(cfg,'Latindex') & cfg.Latindex == true & floor(cfg.clustervars/2) == cfg.clustervars/2
    Latmode = questdlg('Sequence of L and R?','L/R sequence','Cancel','LRLR..','LL..RR..','LL..RR..');
    if strcmp(Latmode,'LL..RR..')
        datatmp = [];
        varstmp = [];
        for i = 1:cfg.numclusters
            startL = (i-1)*cfg.clustervars + 1;
            endL = i * cfg.clustervars - cfg.clustervars/2;
            startR = endL + 1;
            endR = endL + cfg.clustervars/2;
            datatmp = [datatmp ((data(:,startL:endL)-data(:,startR:endR))./(data(:,startL:endL) + data(:,startR:endR)))]; %#ok<AGROW>
            for j = 0:(cfg.clustervars/2-1)
                varstmp = [varstmp  cellstr([header.vars{startL+j} '-' header.vars{startR+j}])]; %#ok<AGROW>
            end
        end
        data = datatmp;
        header.vars = varstmp;
        cfg.clustervars = cfg.clustervars/2;
        cfg.Latindex = 1;
        clearvars varstmp datatmp i j
    elseif strcmp(Latmode,'LRLR..')
        data = (data(:,1:2:end) - data(:,2:2:end)) ./ (data(:,1:2:end) + data(:,2:2:end));
        varstmp = [];
        for j = 1:2:size(header.vars,2)
            varstmp = [varstmp  cellstr([header.vars{j} '-' header.vars{j+1}])]; %#ok<AGROW>
        end
        header.vars = varstmp;
        cfg.clustervars = cfg.clustervars/2;
        cfg.Latindex = 1;
        clearvars varstmp j
    else
        cfg.Latindex = 0;
    end
end

% Calculate z-values
if isfield(cfg,'zvalues') & cfg.zvalues == true
    selection = listdlg('PromptString','Z-value measures:','SelectionMode','multiple', ...
        'InitialValue',(1:size(measures,1)),'ListString',measures,'ListSize',[300 400]);
    for i = selection
        for j = 1:cfg.clustervars
            index = find(~isnan(data(:,(i-1)*cfg.clustervars + j)));
            data(index,(i-1)*cfg.clustervars + j) = (data(index,(i-1)*cfg.clustervars + j) - ...
                mean(data(index,(i-1)*cfg.clustervars + j))) / std(data(index,(i-1)*cfg.clustervars + j));
            header.vars{1,(i-1)*cfg.clustervars + j} = [header.vars{1,(i-1)*cfg.clustervars + j} '_Z'];
        end
        measures{i,1} = [measures{i,1} '_Z']; %#ok<AGROW>
    end
end
clearvars i j selection

% Calculate Log
if isfield(cfg,'Logdata') & cfg.Logdata == true
    selection = listdlg('PromptString','Log measures:','SelectionMode','multiple', ...
        'InitialValue',(1:size(measures,1)),'ListString',measures,'ListSize',[300 400]);
    for i = selection
        for j = 1:cfg.clustervars
            index = find(~isnan(data(:,(i-1)*cfg.clustervars + j)));
            data(index,(i-1)*cfg.clustervars + j) = log(data(index,(i-1)*cfg.clustervars + j));
            header.vars{1,(i-1)*cfg.clustervars + j} = [header.vars{1,(i-1)*cfg.clustervars + j} '_Log'];
        end
        measures{i,1} = [measures{i,1} '_Log']; %#ok<AGROW>
    end
end
clearvars i j selection

% Calculate Logit
if isfield(cfg,'Logitdata') & cfg.Logitdata == true
    selection = listdlg('PromptString','Logit measures:','SelectionMode','multiple', ...
        'InitialValue',(1:size(measures,1)),'ListString',measures,'ListSize',[300 400]);
    for i = selection
        for j = 1:cfg.clustervars
            index = find(~isnan(data(:,(i-1)*cfg.clustervars + j)));
            data(index,(i-1)*cfg.clustervars + j) = lab_logit(data(index,(i-1)*cfg.clustervars + j));
            header.vars{1,(i-1)*cfg.clustervars + j} = [header.vars{1,(i-1)*cfg.clustervars + j} '_Logit'];
        end
        measures{i,1} = [measures{i,1} '_Logit']; %#ok<AGROW>
    end
end
clearvars i j selection

% Select subjects
[selection] = listdlg('PromptString','Select subjects:','SelectionMode','multiple', ...
    'InitialValue',(1:size(header.subjects,1)),'ListString',header.subjects,'ListSize',[300 400]);
if isempty(selection)
    data = [];
    result = [];
    factors = [];
    header = [];
    return
end
header.subjects = header.subjects(selection,:);
data = data(selection,:);
result = result(selection,:);
if ~isempty(factors)
    factors = factors(selection,:);
end
if exist('variables','var')
    header.variables = variables;
end
if exist('measures','var')
    header.measures = measures;
end

% write xls-file
if isfield(cfg,'file') & ~isempty(cfg.file) & isfield(cfg,'path') & ~isempty(cfg.path)
    [~,~,~,filenameS] = lab_filename(cfg.file);
    xlsout{1,1} = ['C' num2str(cfg.clustervars) ' R' num2str(cfg.numresults)];
    xlsout = [xlsout header.vars header.result];
    xlsout = cat(1,xlsout,[header.subjects num2cell(data) num2cell(result)]);
    xlsout = xlsout';
    if size(xlsout,2) > 255
        fileout = fullfile(cfg.path,[filenameS '_statistics.xlsx']);
    else
        fileout = fullfile(cfg.path,[filenameS '_statistics.xls']);
    end
    lab_write_xls(fileout,xlsout);
end