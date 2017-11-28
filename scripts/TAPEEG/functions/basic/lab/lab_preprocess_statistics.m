% Helper function for lab_read_statistics
%
% written by F. Hatz 2012

function [data,header,result,factors,cfg] = lab_preprocess_statistics(inputdata,cfg,doresults,dofactors,doselection,dowriting)

data = [];
result = [];
factors = [];
header = [];

if ~exist('doselection','var') | isempty(doselection)
    doselection = 1;
end
if ~exist('dofactors','var') | isempty(dofactors)
    dofactors = 1;
end
if ~exist('doresults','var') | isempty(doresults)
    doresults = 1;
end

if ~exist('inputdata','var') | isempty(inputdata)
    disp('Input file is invalid')
    return
else
    inputdata = inputdata';
end

if ~exist('cfg','var') | ~isfield(cfg,'clustervars')
    cfg.clustervars = [];
    cfg.numresults = [];
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
end
if ~isfield(cfg,'numresults') | isempty(cfg.numresults)
    cfg.numresults = 0;
end

if cfg.numresults > 0
    cfg.ResultNames = {};
    Result = inputdata(2:end,end-cfg.numresults+1:end);
    for i = 1:cfg.numresults
        if ~isnumeric(Result{1,i})
            [tmp,~,tmp2] = unique(Result(:,i),'stable');
            Result(:,i) = num2cell(tmp2(:));
            cfg.ResultNames{i} = tmp(:);
        else
            tmp = Result(:,i);
            for j = 1:length(tmp)
                if isnumeric(tmp{j})
                    tmp{j} = num2str(tmp{j});
                end
            end
            tmp = unique(tmp,'stable');
            cfg.ResultNames{i} = tmp(:);
        end
        
    end
    inputdata(2:end,end-cfg.numresults+1:end) = Result;
end

try
    data = cell2mat(inputdata(2:end,2:end));
catch %#ok<CTCH>
    disp('Input file has invalid / non numeric data')
    return
end

if isempty(cfg.clustervars)
    Prompt = {};
    Formats = {};
    Prompt(end+1,:) = {'Number of variables per measure', 'clustervars'};
    Formats(end+1,1).type = 'edit';
    Formats(end,1).format = 'integer';
    Formats(end,1).limits = [0 9999];
    Formats(end,1).size = 60;
    [cfg,Cancelled] = inputsdlg(Prompt,'Variables per measure',Formats,cfg);
    if Cancelled == 1
        data = [];
        return
    else
        cfg.numresults = cfg.numresults + mod((size(data,2)-cfg.numresults),cfg.clustervars);
    end
end
if ~isfield(cfg,'clustervars') | isempty(cfg.clustervars) | size(data,2) < cfg.clustervars
    data = [];
    disp('Abort: Number of variables per measure and data-file not matching')
    return
end

if cfg.clustervars > 1
    if isfield(cfg,'numresults') & ~isempty(cfg.numresults)
        cfg.numclusters = floor((size(data,2)-cfg.numresults)/cfg.clustervars);
    else
        cfg.numclusters = floor(size(data,2)/cfg.clustervars);
    end
    cfg.numresults = size(data,2) - (cfg.numclusters*cfg.clustervars);
else
    cfg.clustervars = 1;
    cfg.numclusters = size(data,2) - cfg.numresults;
end
if ~isfield(cfg,'clustervars2')
    cfg.clustervars2 = 1;
elseif mod(cfg.numclusters,cfg.clustervars2) ~= 0
    cfg.clustervars2 = 1;
end
cfg.numclusters = cfg.numclusters / cfg.clustervars2;

% find possible results
if ~isfield(cfg,'numresults') | isempty(cfg.numresults)
    cfg.numresults = 0;
end
if cfg.numresults > 0
    resultnames = inputdata(1,end-cfg.numresults+1:end);
    result = data(:,end-cfg.numresults+1:end);
    varnames = inputdata(1,2:end-cfg.numresults);
    data = data(:,1:end-cfg.numresults);
elseif doresults == 1
    flagresultall = 1; %#ok<NASGU>
    resultnames = inputdata(1,2:end);
    result = data;
    varnames = inputdata(1,2:end);
else
    resultnames = {};
    result = [];
    varnames = inputdata(1,2:end);
end
if ~isempty(resultnames)
    resultnames = regexprep(resultnames,{'::',':',';','/','\',','},'_');
    cfg.selectresults = length(resultnames);
else
    cfg.selectresults = [];
end
cfg.selectfactors = [];

% find names of measures and variables
if cfg.clustervars > 1
    for i = 1:cfg.numclusters*cfg.clustervars2
        tmp = strfind(varnames{1,(i-1)*cfg.clustervars+1},'_');
        if cfg.clustervars < 2
            measures{1,i} = varnames{1,(i-1)*cfg.clustervars+1}; %#ok<AGROW>
        elseif ~isempty(tmp)
            measures{1,i} = varnames{1,(i-1)*cfg.clustervars+1}(1:tmp(end)-1); %#ok<AGROW>
        else
            measures{1,i} = ['Measure' num2str(i)]; %#ok<AGROW>
        end
    end
    clearvars tmp i
    for i = 1:cfg.clustervars
        tmp = strfind(varnames{1,i},'_');
        if ~isempty(tmp)
            variables{1,i} = varnames{1,i}(tmp(end)+1:end); %#ok<AGROW>
        else
            variables{1,i} = ['Var' num2str(i)]; %#ok<AGROW>
        end
    end
    clearvars tmp i
    if cfg.clustervars2 > 1
        varnames2 = measures;
        clearvars measures
        for i = 1:cfg.numclusters
            tmp = strfind(varnames2{1,(i-1)*cfg.clustervars2+1},'_');
            if cfg.clustervars2 < 2
                measures{1,i} = varnames2{1,(i-1)*cfg.clustervars2+1}; %#ok<AGROW>
            elseif ~isempty(tmp)
                measures{1,i} = varnames2{1,(i-1)*cfg.clustervars2+1}(1:tmp(end)-1); %#ok<AGROW>
            else
                measures{1,i} = ['Measure' num2str(i)]; %#ok<AGROW>
            end
        end
        clearvars tmp i
        for i = 1:cfg.clustervars2
            tmp = strfind(varnames2{1,i},'_');
            if ~isempty(tmp)
                variables2{1,i} = varnames2{1,i}(tmp(end)+1:end); %#ok<AGROW>
            else
                variables2{1,i} = ['Var' num2str(i)]; %#ok<AGROW>
            end
        end
        clearvars tmp i
        variables2 = regexprep(variables2,{'::',':',';','/','\',','},'_');
        clearvars varnames2
    else
        variables2 = {};
    end
    measures = regexprep(measures,{'::',':',';','/','\',','},'_');
    variables = regexprep(variables,{'::',':',';','/','\',','},'_');
else
    measures = inputdata(1,2:end-cfg.numresults)';
    measures = regexprep(measures,{'::',':',';','/','\',','},'_');
    variables = {};
    variables2 = {};
end
if doselection == 2
    cfg.selectmeasures = 1;
else
    cfg.selectmeasures = 1:length(measures);
end
if ~isempty(variables)
    cfg.selectvariables = 1:length(variables);
else
    cfg.selectvariables = 1;
end
if ~isempty(variables2)
    cfg.selectvariables2 = 1:length(variables2);
else
    cfg.selectvariables2 = 1;
end

% set names of subject / trials
subjects = inputdata(2:end,1);
subjects = regexprep(subjects,' ','');
cfg.selectsubjects = 1:length(subjects);

% rearrange data
data = reshape(data,[size(data,1) cfg.clustervars cfg.clustervars2 cfg.numclusters]);
    
% Ask for settings of input data
[cfg,skipprocessing] = lab_set_preprocess_statistics(cfg,measures,variables,variables2,resultnames,subjects,doselection,dofactors,dowriting);
if skipprocessing == 1
    data = [];
    return
else
    pause(0.2);
end

% correct selections
cfg.selectfactors = setdiff(cfg.selectfactors,1);
if ~isempty(cfg.selectfactors)
    cfg.selectfactors = cfg.selectfactors - 1;
    cfg.selectfactors = setdiff(cfg.selectfactors,cfg.selectresults);
end
if cfg.clustervars == 1 & exist('flagresultall','var') & doresults == 1
    cfg.selectmeasures = setdiff(cfg.selectmeasures,cfg.selectresults);
    cfg.selectmeasures = setdiff(cfg.selectmeasures,cfg.selectfactors);
end

% Reduce data to selections
data = data(cfg.selectsubjects,cfg.selectvariables,cfg.selectvariables2,cfg.selectmeasures);
header.vars = {};
for i = 1:length(cfg.selectmeasures)
    for j = 1:length(cfg.selectvariables2)
        for m = 1:length(cfg.selectvariables)
            if ~isempty(variables) & ~isempty(variables2)
                header.vars{m,j,i} = [measures{cfg.selectmeasures(i)} '_' ...
                    variables2{cfg.selectvariables2(j)} '_' variables{cfg.selectvariables(m)}];
            elseif ~isempty(variables)
                header.vars{m,j,i} = [measures{cfg.selectmeasures(i)} '_' variables{cfg.selectvariables(m)}];
            else
                header.vars{m,j,i} = measures{cfg.selectmeasures(i)};
            end
        end
    end
end
measures = measures(cfg.selectmeasures);
cfg.numclusters = length(cfg.selectmeasures);
if ~isempty(variables2)
    variables2 = variables2(cfg.selectvariables2);
    cfg.clustervars2 = length(cfg.selectvariables2);
end
if ~isempty(variables)
    variables = variables(cfg.selectvariables);
    cfg.clustervars = length(cfg.selectvariables);
end
subjects = subjects(cfg.selectsubjects);
header.subjects = subjects(:);
if ~isempty(cfg.selectfactors);
    factors = result(cfg.selectsubjects,cfg.selectfactors);
    header.factors = resultnames(cfg.selectfactors);
else
    header.factors = {};
end
if ~isempty(cfg.selectresults) & ~isempty(result)
    result = result(cfg.selectsubjects,cfg.selectresults);
    cfg.numresults = length(cfg.selectresults);
    header.result = resultnames(cfg.selectresults);
    if isfield(cfg,'ResultNames')
        cfg.ResultNames = cfg.ResultNames(cfg.selectresults);
    end
else
    header.result = {};
end

% Calculate z-values
if isfield(cfg,'zvalues') & cfg.zvalues == true
    selection = listdlg('PromptString','Z-value measures:','SelectionMode','multiple', ...
        'InitialValue',(1:length(measures)),'ListString',measures,'ListSize',[300 400]);
    for i = selection
        for j = 1:cfg.clustervars2
            for m = 1:cfg.clustervars
                index = find(~isnan(data(:,m,j,i)));
                data(index,m,j,i) = (data(index,m,j,i) - ...
                    mean(data(index,m,j,i))) / std(data(index,m,j,i));
                header.vars{m,j,i} = ['Z_' header.vars{m,j,i}];
            end
        end
        measures{i} = ['Z_' measures{i}];
    end
end
clearvars i j m selection

% Calculate Log
if isfield(cfg,'Logdata') & cfg.Logdata == true
    selection = listdlg('PromptString','Log measures:','SelectionMode','multiple', ...
        'InitialValue',(1:length(measures)),'ListString',measures,'ListSize',[300 400]);
    for i = selection
        for j = 1:cfg.clustervars2
            for m = 1:cfg.clustervars
                index = find(~isnan(data(:,m,j,i)));
                data(index,m,j,i) = log(data(index,m,j,i));
                header.vars{m,j,i} = ['Log_' header.vars{m,j,i}];
            end
        end
        measures{i} = ['Log_' measures{i}];
    end
end
clearvars i j m selection

% Calculate Logit
if isfield(cfg,'Logitdata') & cfg.Logitdata == true
    selection = listdlg('PromptString','Logit measures:','SelectionMode','multiple', ...
        'InitialValue',(1:length(measures)),'ListString',measures,'ListSize',[300 400]);
    for i = selection
        for j = 1:cfg.clustervars2
            for m = 1:cfg.clustervars
                index = find(~isnan(data(:,m,j,i)));
                data(index,m,j,i) = lab_logit(data(index,m,j,i));
                header.vars{m,j,i} = ['Logit_' header.vars{m,j,i}];
            end
        end
        measures{i} = ['Logit_' measures{i}];
    end
end
clearvars i j selection

% Calculate ratio
if isfield(cfg,'calcratio') & cfg.calcratio == true
    selection = listdlg('PromptString','Ratio first measure:','SelectionMode','single', ...
        'ListString',measures,'ListSize',[300 400]);
    selection2 = listdlg('PromptString','Ratio second measure:','SelectionMode','single', ...
        'ListString',measures,'ListSize',[300 400]);
    if ~isempty(selection) & ~isempty(selection2)
        data = real(log(data(:,:,:,selection) ./ data(:,:,:,selection2)));
        measures = {['LogRatio_' measures{selection} '_' measures{selection2}]};
        for i = 1:cfg.clustervars2
            header.vars = {};
            for j = 1:cfg.clustervars
                if ~isempty(variables) & ~isempty(variables2)
                    header.vars{j,i} = [measures{1} '_' variables2{i} '_' variables{j}];
                elseif ~isempty(variables)
                    header.vars{j,i} = [measures{1} '_' variables{j}];
                else
                    header.vars{j,i} = measures{1};
                end
            end
        end
        cfg.numclusters = 1;
    end
end
clearvars i j selection selection2

% Calculate lateralization index
if isfield(cfg,'Latindex') & cfg.Latindex == true & floor(cfg.clustervars/2) == cfg.clustervars/2
    Latmode = questdlg('Sequence of L and R?','L/R sequence','Cancel','LRLR..','LL..RR..','LL..RR..');
    if strcmp(Latmode,'LL..RR..')
        Left = 1:cfg.clustervars/2;
        Right = cfg.clustervars/2+1:cfg.clustervars;
        cfg.Latindex = 1;
    elseif strcmp(Latmode,'LRLR..')
        Left = 1:2:cfg.clustervars;
        Right = 2:2:cfg.clustervars;
        cfg.Latindex = 1;
    else
        cfg.Latindex = 0;
    end
    if cfg.Latindex == 1
        data = (data(:,Left,:,:) - data(:,Right,:,:))./(data(:,Left,:,:) + data(:,Right,:,:));
        header.vars = {};
        for i = 1:length(measures)
            for j = 1:cfg.clustervars2
                for m = 1:cfg.clustervars/2
                    if ~isempty(variables) & ~isempty(variables2)
                        header.vars{m,j,i} = [measures{i} '_' variables2{j} '_' variables{Left(m)} '-' variables{Right(m)}];
                    elseif ~isempty(variables)
                        header.vars{m,j,i} = [measures{i} '_' variables{Left(m)} '-' variables{Right(m)}];
                    end
                end
            end
        end
        cfg.clustervars = cfg.clustervars/2;
        clearvars  i j m
    end
end

% Calculate factors
if isfield(cfg,'calcfactors') & cfg.calcfactors == true
    if ~isfield(cfg,'fcomb')
        cfg.factorscalc = {};
        if length(measures) > 2
            cfg.factorscalc{1,1} = 'V001 + V002 + V003';
        elseif length(measures) > 1
            cfg.factorscalc{1,1} = 'V001 + V002';
        else
            cfg.factorscalc{1,1} = 'V001 / 3';
        end
    end
    cfg.listmeasure = {};
    for i = 1:length(measures)
        cfg.listmeasure{i,1} = ['V' num2str(i,'%03d')];
    end
    cfg.listmeasure = [cfg.listmeasure measures(:)'];
    
    Prompt = cell(0,2);
    Formats = [];
    
    Prompt(end+1,:) = {'factors','factorscalc'};
    Formats(end+1,1).type = 'table';
    Formats(end,1).format = 'table';
    Formats(end,1).items = {{'Formula'},[],{210}};
    Formats(end,1).callback = {@lab_table_dialog,[],'factorscalc',{'Formula'}};
    Formats(end,1).size = [250 150];
    Formats(end,1).enable = 'inactive';
    
    Prompt(end+1,:) = {'measures','listmeasure'};
    Formats(end+1,1).type = 'table';
    Formats(end,1).format = 'table';
    Formats(end,1).items = {{'string','measure'},[],{50 125}};
    Formats(end,1).size = [250 150];
    Formats(end,1).enable = 'inactive';
    
    [cfg,Cancelled] = inputsdlg(Prompt,'Create factors',Formats,cfg,2);
    if Cancelled == 1
        data = [];
        return
    else
        cfg = rmfield(cfg,'listmeasure');
    end
    
    vars = {};
    datanew = [];
    for nfact = 1:length(cfg.factorscalc)
        measuresNew{nfact,1} = ['Factor' num2str(nfact,'%03d')]; %#ok<AGROW>
        for j = 1:cfg.clustervars2
            for i = 1:cfg.clustervars
                formula = cfg.factorscalc{nfact};
                for m = 1:length(measures)
                    formula = regexprep(formula,['V' num2str(m,'%03d')],['data(:,' num2str(i) ',' num2str(j) ',' num2str(m) ')']);
                end
                formula = regexprep(formula,'/','./');
                formula = regexprep(formula,'*','.*');
                tmp = eval(formula);
                if size(tmp,1) == size(data,1)
                    datanew(:,i,j,nfact) = tmp; %#ok<AGROW>
                elseif isnumeric(tmp) & length(tmp) == 1
                    datanew(:,i,j,nfact) = repmat(tmp,size(data,1),1); %#ok<AGROW>
                else
                    datanew(:,i,j,nfact) = NaN(size(data,1),1); %#ok<AGROW>
                end
                if ~isempty(variables) & ~isempty(variables2)
                    vars{i,j,nfact} = [measuresNew{nfact,1} '_' variables2{j} '_' variables{i}]; %#ok<AGROW>
                elseif ~isempty(variables)
                    vars{i,j,nfact} = [measuresNew{nfact,1} '_' variables{i}]; %#ok<AGROW>
                else
                    vars{i,j,nfact} = measuresNew{nfact,1}; %#ok<AGROW>
                end
            end
        end
    end
    header.vars = vars;
    data = datanew;
    measures = measuresNew;
    cfg.numclusters = length(measures);
    clearvars vars datanew
end

if exist('measures','var')
    header.measures = measures(:)';
end
if ~isempty(variables)
    header.variables = variables(:)';
end
if ~isempty(variables2)
    header.variables2 = variables2(:)';
end

% Exclude trials with invalid entries
if isfield(cfg,'excludeNaN') & cfg.excludeNaN == true
    tmp = find(max(isnan(reshape(data,size(data,1),size(data,2)*size(data,3)*size(data,4))),[],2));
    if ~isempty(result)
        tmp = union(tmp,find(max(isnan(result),[],2)));
    end
    if ~isempty(factors)
        tmp = union(tmp,find(max(isnan(factors),[],2)));
    end
    tmp = setdiff(1:size(data,1),tmp);
    data = data(tmp,:,:,:);
    if ~isempty(result)
        result = result(tmp,:);
    end
    if ~isempty(factors)
        factors = factors(tmp,:);
    end
    header.subjects = header.subjects(tmp,1);
    clearvars tmp
end

% Change order of measures / variables
if isfield(cfg,'changeorder') & ~isempty(cfg.changeorder)
    VARS = {subjects,variables,variables2,measures};
    Map = [1 4 3 2];
    changeorder = [Map(cfg.changeorder.subjects) Map(cfg.changeorder.variables) Map(cfg.changeorder.variables2) Map(cfg.changeorder.measures)];
    variables = VARS{changeorder(2)}(:)';
    variables2 = VARS{changeorder(3)}(:)';
    measures = VARS{changeorder(4)}(:)';
    subjects = VARS{changeorder(1)}(:);
    data = permute(data,changeorder);
    cfg.clustervars2 = size(data,3);
    cfg.clustervars = size(data,2);
    header.vars = {};
    for i = 1:length(measures)
        for j = 1:cfg.clustervars2
            for m = 1:cfg.clustervars
                if ~isempty(variables) & ~isempty(variables2)
                    header.vars{m,j,i} = [measures{i} '_' variables2{j} '_' variables{m}];
                elseif ~isempty(variables)
                    header.vars{m,j,i} = [measures{i} '_' variables{m}];
                end
            end
        end
    end
        
    variables2 = variables2(:);
    Skip = false;
    for i = 1:size(variables2,1)
        tmp = strfind(variables2{i,1},'_');
        if ~isempty(tmp)
            Vars{i,1} = variables2{i,1}(1:tmp(end)-1); %#ok<AGROW>
            Vars2{i,1} = variables2{i,1}(tmp(end)+1:end); %#ok<AGROW>
        else
            Skip = true;
        end
    end
    if Skip == false & ~isempty(variables2)
        [~,m,~]=unique(Vars,'last');
        m = sort(m);
        if length(m) == 1 | (m(1) > 1 & (m(2)-m(1)) == m(1))
            cfg.clustervars2 = m(1);
            variables2 = Vars2(1:cfg.clustervars2);
            Vars = Vars(1:cfg.clustervars2:end);
            tmp = {};
            for i = 1:length(measures)
                for j = 1:length(Vars)
                    tmp{end+1,1} = [measures{i} '_' Vars{j}]; %#ok<AGROW>
                end
            end
            measures = tmp;
        end
        clearvars tmp m i
    end
    
    header.variables = variables(:)';
    header.variables2 = variables2(:)';
    header.measures = measures(:)';
    header.subjects = subjects(:);
end

% Combine Measures/Variables and Subjects/Trials
if isfield(cfg,'combinevars') & isfield(cfg.combinevars,'Variable') & ~isempty(cfg.combinevars.Variable)
    if cfg.combinevars.Variable == 1
        if ~isempty(variables2)
            data = permute(data,[4 1 2 3]);
            data = reshape(data,[size(data,1)*size(data,2) size(data,3) 1 size(data,4)]);
            newvars = measures;
            measures = variables2;
            variables2 = {};
            cfg.clustervars2 = 1;
        else
            data = permute(data,[4 1 3 2]);
            data = reshape(data,[size(data,1)*size(data,2) size(data,3) 1 size(data,4)]);
            newvars = measures;
            measures = variables;
            variables = {};
            variables2 = {};
            cfg.clustervars = 1;
            cfg.clustervars2 = 1;
        end
    elseif cfg.combinevars.Variable == 2
        if ~isempty(variables2)
            data = permute(data,[3 1 2 4]);
            data = reshape(data,[size(data,1)*size(data,2) size(data,3) 1 size(data,4)]);
            newvars = variables2;
            variables2 = {};
            cfg.clustervars2 = 1;
        else
            data = permute(data,[2 1 3 4]);
            data = reshape(data,[size(data,1)*size(data,2) size(data,3) 1 size(data,4)]);
            newvars = variables;
            variables = {};
            variables2 = {};
            cfg.clustervars = 1;
            cfg.clustervars2 = 1;
        end
    elseif cfg.combinevars.Variable == 3
        if ~isempty(variables2)
            data = permute(data,[2 1 3 4]);
            data = reshape(data,[size(data,1)*size(data,2) size(data,3) 1 size(data,4)]);
            newvars = variables;
            variables = variables2;
            variables2 = {};
            cfg.clustervars = cfg.clustervars2;
            cfg.clustervars2 = 1;
        else
            newvars = {};
        end
    end
    subjectstmp = cell(0,1);
    for j = 1:length(subjects)
        for i = 1:length(newvars)
            subjectstmp{end+1,1} = [subjects{j} '_' newvars{i}]; %#ok<AGROW>
        end
    end
    if ~isempty(result)
        resulttmp = repmat(result,[1 1 length(newvars)]);
        resulttmp = permute(resulttmp,[3 1 2]);
        result = reshape(resulttmp,[size(resulttmp,1)*size(resulttmp,2) size(resulttmp,3)]);
        clearvars resulttmp
    end
    if ~isempty(factors)
        factorstmp = repmat(factors,[1 1 length(newvars)]);
        factorstmp = permute(factorstmp,[3 1 2]);
        factors = reshape(factorstmp,[size(factorstmp,1)*size(factorstmp,2) size(factorstmp,3)]);
        clearvars factorstmp
    end
    header.vars = {};
    for i = 1:length(measures)
        for j = 1:cfg.clustervars2
            for m = 1:cfg.clustervars
                if ~isempty(variables) & ~isempty(variables2)
                    header.vars{m,j,i} = [measures{i} '_' variables2{j} '_' variables{m}];
                elseif ~isempty(variables)
                    header.vars{m,j,i} = [measures{i} '_' variables{m}];
                else
                    header.vars{m,j,i} = measures{i};
                end
            end
        end
    end
    header.subjects = subjectstmp;
    if isfield(header,'variables2')
        header = rmfield(header,'variables2');
    end
    clearvars subjectstmp i j m
else
    cfg.numclusters = cfg.numclusters*cfg.clustervars2;
end

% rearrange order of outcomes
if isfield(cfg,'reorderresults') & cfg.reorderresults == true & ~isempty(cfg.selectresults)
    for i = 1:length(cfg.selectresults)
        Names = cfg.ResultNames{cfg.selectresults(i)};
        Order = lab_rearrange_names(Names,['Re-arrange ' header.result{cfg.selectresults(i)}]);
        Result = -result(:,cfg.selectresults(i));
        for j = 1:length(Order)
            Result(Result==-j) = Order(j);
        end
        result(:,cfg.selectresults(i)) = Result;
        cfg.ResultNames{cfg.selectresults(i)} = cfg.ResultNames{cfg.selectresults(i)}(Order);
    end
end

data = reshape(data,[size(data,1) size(data,2)*size(data,3)*size(data,4)]);
header.vars = header.vars(:)';