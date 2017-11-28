function [cfg,skipprocessing] = lab_set_preprocess_statistics(cfg,measures,variables,variables2,resultnames,subjects,doselection,dofactors,dowriting)

skipprocessing = 0;

disp('Settings for data')
if ~isfield(cfg,'calcratio')
    cfg.calcratio = false;
    cfg.excludeM = false;
    cfg.excludeV = false;
    cfg.Latindex = false;
    cfg.zvalues = false;
    cfg.Logdata = false;
    cfg.Logitdata = false;
    cfg.excludeNaN = false;
    cfg.calcfactors = false;
    cfg.combinevars = [];
    cfg.changeorder = [];
end

if isfield(cfg,'selectsubjects') & ~isempty(cfg.selectsubjects);
    cfg.subjects = subjects(cfg.selectsubjects);
end
if dowriting == 1
    if isfield(cfg,'file') & isfield(cfg,'path')
        if isfield(cfg,'subfolder') & ~isempty(cfg.subfolder)
            if ~exist(fullfile(cfg.path,cfg.subfolder),'dir')
                mkdir(fullfile(cfg.path,cfg.subfolder));
            end
            cfg.filename = fullfile(fullfile(cfg.path,cfg.subfolder),[cfg.file '.xlsx']);
        else
            cfg.filename = fullfile(cfg.path,[cfg.file '.xlsx']);
        end
    else
        cfg.filename = '';
    end
end
if ~isfield(cfg,'file')
    cfg.file = '';
end
if ~isfield(cfg,'path')
    cfg.path = pwd;
end

Prompt = cell(0,2);
Formats = {};

if dowriting == 1
    Prompt(end+1,:) = {'File to store','filename'};
    Formats(end+1,1).type = 'edit';
    Formats(end,1).format = 'file';
    Formats(end,1).items = {'*.xlsx;*.xls','Excel File (*.xlsx)'};
    Formats(end,1).limits = [1 0];
    Formats(end,1).size = [-1 0];
    Formats(end,1).span = [1 6];
else
    Prompt(end+1,:) = {'Filename','file'};
    Formats(end+1,1).type = 'edit';
    Formats(end,1).format = 'text';
    Formats(end,1).size = 200;
    Formats(end,1).span = [1 6];
end

Prompt(end+1,:) = {'Measures','selectmeasures'};
Formats(end+1,1).type = 'list';
if doselection == 2
    Formats(end,1).style = 'popupmenu';
    Formats(end,1).items = measures;
    Formats(end,1).size = 140;
else
    Formats(end,1).style = 'listbox';
    Formats(end,1).items = measures;
    Formats(end,1).limits = [0 2]; % multi-select
    Formats(end,1).size = [140 160];
    Formats(end,1).span = [1 2];
end

if ~isempty(variables) & ~isempty(variables2)
    Prompt(end+1,:) = {'1. Variables','selectvariables2'};
    Formats(end+1,1).type = 'list';
    Formats(end,1).style = 'listbox';
    Formats(end,1).items = variables2;
    Formats(end,1).limits = [0 2]; % multi-select
    Formats(end,1).size = [140 160];
    Formats(end,1).span = [1 2];
    
    Prompt(end+1,:) = {'2. Variables','selectvariables'};
    Formats(end+1,1).type = 'list';
    Formats(end,1).style = 'listbox';
    Formats(end,1).items = variables;
    Formats(end,1).limits = [0 2]; % multi-select
    Formats(end,1).size = [140 160];
    Formats(end,1).span = [1 2];
elseif  ~isempty(variables)
    Prompt(end+1,:) = {'Variables','selectvariables'};
    Formats(end+1,1).type = 'list';
    Formats(end,1).style = 'listbox';
    Formats(end,1).items = variables;
    Formats(end,1).limits = [0 2]; % multi-select
    Formats(end,1).size = [140 160];
    Formats(end,1).span = [1 2];
    
    Formats(end+1,1).type = 'none';
    Formats(end,1).span = [1 2];
else
    Formats(end+1,1).type = 'none';
    Formats(end,1).span = [1 4];
end

Formats(end+1,1).type = 'none';
Formats(end,1).span = [1 6];

if ~isempty(resultnames)
    Prompt(end+1,:) = {'Outcome ','selectresults'};
    Formats(end+1,1).type = 'list';
    Formats(end,1).style = 'listbox';
    Formats(end,1).items = resultnames;
    Formats(end,1).limits = [0 2]; % multi-select
    Formats(end,1).size = [140 120];
    Formats(end,1).span = [1 2];
    
    if dofactors == 1
        Prompt(end+1,:) = {'Factors  ','selectfactors'};
        Formats(end+1,1).type = 'list';
        Formats(end,1).style = 'listbox';
        Formats(end,1).items = [cellstr('none') resultnames];
        Formats(end,1).limits = [0 2]; % multi-select
        Formats(end,1).size = [140 120];
        Formats(end,1).span = [1 2];
    else
        Formats(end+1,1).type = 'none';
        Formats(end,1).span = [1 2];
    end
    
    Formats(end+1,1).type = 'none';
    Formats(end,1).span = [1 2];
end

Formats(end+1,1).type = 'none';
Formats(end,1).span = [1 6];

Prompt(end+1,:) = {'Z-values' 'zvalues'};
Formats(end+1,1).type = 'check';

Prompt(end+1,:) = {'Log data' 'Logdata'};
Formats(end+1,1).type = 'check';

Prompt(end+1,:) = {'Logit data' 'Logitdata'};
Formats(end+1,1).type = 'check';

Formats(end+1,1).type = 'none';
Formats(end,1).span = [1 3];

Prompt(end+1,:) = {'Calculate ratio' 'calcratio'};
Formats(end+1,1).type = 'check';
Formats(end,1).size = [-1 -1];

Prompt(end+1,:) = {'Calculate factors' 'calcfactors'};
Formats(end+1,1).type = 'check';

if cfg.clustervars > 1
    Prompt(end+1,:) = {'Lateralization index' 'Latindex'};
    Formats(end+1,1).type = 'check';
end

if ~isempty(variables)
    Prompt(end+1,:) = {'Change order of variables' 'changeorder'};
    Formats(end+1,1).type = 'check';
    Formats(end,1).callback = {@change_order,'changeorder','changeorder',variables,variables2};
    Formats(end,1).span = [1 6];
end

Prompt(end+1,:) = {'Combine variables and subjects/trials' 'combinevars'};
Formats(end+1,1).type = 'check';
Formats(end,1).callback = {@combine_vars,'@ALL','@ALL',variables,variables2};
Formats(end,1).span = [1 6];

if ~isempty(resultnames)
    Prompt(end+1,:) = {'Re-order outcomes' 'reorderresults'};
    Formats(end+1,1).type = 'check';
    Formats(end,1).span = [1 6];
end

Formats(end+1,1).type = 'none';
Formats(end,1).span = [1 6];

Prompt(end+1,:) = {'Subjects/Trials','subjects'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'vector';
Formats(end,1).items = subjects;
Formats(end,1).limits = [-1 1];
Formats(end,1).size = [370 20];
Formats(end,1).span = [1 6];

Prompt(end+1,:) = {'Import subjects',''};
Formats(end+1,1).type = 'button';
Formats(end,1).style = 'pushbutton';
Formats(end,1).size = [100 25];
Formats(end,1).callback = {@load_subjects,'subjects','subjects',subjects};

Prompt(end+1,:) = {'Save subjects',''};
Formats(end+1,1).type = 'button';
Formats(end,1).style = 'pushbutton';
Formats(end,1).size = [100 25];
Formats(end,1).callback = {@save_subjects,[],'subjects'};

Prompt(end+1,:) = {'Exclude trials with invalid entries' 'excludeNaN'};
Formats(end+1,1).type = 'check';
Formats(end,1).span = [1 4];

[cfg,Cancelled] = inputsdlg(Prompt,'Select data statistics',Formats,cfg,3);
if Cancelled == 1
    skipprocessing = 1;
else
    if verLessThan('matlab','8')
        [~,cfg.selectsubjects] = intersect(subjects,cfg.subjects);
    else
        [~,cfg.selectsubjects] = intersect(subjects,cfg.subjects,'stable');
    end
    if dowriting == 1
        [~,tmp2,~,tmp1] = lab_filename(cfg.filename);
        cfg.file = tmp1;
        cfg.path = tmp2;
        cd(cfg.path);
        clearvars tmp1 tmp2
    end
end

end

function subjectsselect = load_subjects(subjectsselect,subjects)
   [Subjects_file,Subjects_filepath] = uigetfile('.xls;..xlsx','Select file with subjects');
   if isnumeric(Subjects_file)
       return
   end
   SubjectsFile = lab_read_xls(fullfile(Subjects_filepath,Subjects_file));
   if isempty(SubjectsFile)
       return
   end
   if size(SubjectsFile,2) == 1
       SubjectsFile = SubjectsFile';
   end
   if verLessThan('matlab','8')
       [~,selection] = intersect(subjects,SubjectsFile(1,:));
       if isempty(selection)
           [~,selection] = intersect(subjects,SubjectsFile(:,1));
       end
   else
       [~,selection] = intersect(subjects,SubjectsFile(1,:),'stable');
       if isempty(selection)
           [~,selection] = intersect(subjects,SubjectsFile(:,1),'stable');
       end
   end
   if ~isempty(selection)
       subjectsselect = subjects(selection);
   end
end

function save_subjects(subjectsselect)
   [Subjects_file,Subjects_filepath] = uiputfile('subjects.xlsx','Select file to store subjects');
   lab_write_xls(fullfile(Subjects_filepath,Subjects_file),subjectsselect);
end

function changeorder = change_order(changeorder,variables,variables2)
    List = {'Subjects','Measures'};
    if ~isfield(changeorder,'First')
        changeorder.subjects = 1;
        changeorder.measures = 2;
        if ~isempty(variables)
            if ~isempty(variables2)
                changeorder.variables2 = 3;
                changeorder.variables = 4;
            else
                changeorder.variables = 3;
            end
        end
    end
    if ~isempty(variables)
        List = cat(2,List,cellstr('1. Variables'));
        if ~isempty(variables2)
            List = cat(2,List,cellstr('2. Variables'));
        end
    end
    
    Prompt = {};
    Formats = {};
    
    Prompt(end+1,:) = {'Subjects','subjects'};
    Formats(end+1,1).type = 'list';
    Formats(end,1).style = 'popupmenu';
    Formats(end,1).items = List;
    Formats(end,1).callback = {@switch_order,'@ALL','@ALL',1};
    
    Prompt(end+1,:) = {'Measures','measures'};
    Formats(end+1,1).type = 'list';
    Formats(end,1).style = 'popupmenu';
    Formats(end,1).items = List;
    Formats(end,1).callback = {@switch_order,'@ALL','@ALL',2};
    
    if length(List) == 4
        Prompt(end+1,:) = {'1. Variables','variables2'};
        Formats(end+1,1).type = 'list';
        Formats(end,1).style = 'popupmenu';
        Formats(end,1).items = List;
        Formats(end,1).callback = {@switch_order,'@ALL','@ALL',3};
        
        Prompt(end+1,:) = {'2. Variables','variables'};
        Formats(end+1,1).type = 'list';
        Formats(end,1).style = 'popupmenu';
        Formats(end,1).items = List;
        Formats(end,1).callback = {@switch_order,'@ALL','@ALL',4};
    elseif length(List) == 3
        Prompt(end+1,:) = {'1. Variables','variables'};
        Formats(end+1,1).type = 'list';
        Formats(end,1).style = 'popupmenu';
        Formats(end,1).items = List;
        if changeorder.subjects == 4
            changeorder.subjects = 3;
        elseif changeorder.measures == 4
            changeorder.measures = 3;
        elseif changeorder.variables == 4
            changeorder.variables = 3;
        end
        Formats(end,1).callback = {@switch_order,'@ALL','@ALL',3};
    end
    
    [changeorder,Cancelled] = inputsdlg(Prompt,'Change order',Formats,changeorder);
    if Cancelled == 1
        changeorder = [];
    elseif length(List) == 3
        if changeorder.subjects == 3
            changeorder.subjects = 4;
        elseif changeorder.measures == 3
            changeorder.measures = 4;
        elseif changeorder.variables == 3
            changeorder.variables = 4;
        end
        changeorder.variables2 = 3;
    elseif length(List) == 2
        changeorder.variables = 3;
        changeorder.variables2 = 4;
    end
    
    function order = switch_order(order,Nr)
        if isfield(order,'variables2')
            tmp = {'subjects','measures','variables2','variables'};
        elseif isfield(order,'variables')
            tmp = {'subjects','measures','variables'};
        else
            tmp = {'subjects','measures'};
        end
        actorder = order.(tmp{Nr});
        tmp2 = 1:length(tmp);
        for i = 1:length(tmp);
            tmp2(i) = order.(tmp{i});
        end
        misorder = setdiff(1:length(tmp),tmp2);
        tmp2(Nr) = 0;
        tmp2 = find(tmp2 == actorder);
        if ~isempty(tmp2)
            order.(tmp{tmp2}) = misorder;
        end
    end
end

function settings = combine_vars(settings,variables,variables2)
    List = cellstr('Measures');
    if ~isfield(settings,'combinevar') | ~isfield(settings.combinevars,'Variable')
        settings.combinevars.Variable = 1;
    end
    if ~isempty(variables)
        List = cat(2,List,cellstr('1. Variables'));
    end
    if ~isempty(variables2)
        List = cat(2,List,cellstr('2. Variables'));
    end
    [settings.combinevars,Cancelled] = inputsdlg(settings.combinevars,'Select Variable',{List});
    if Cancelled == 1
        settings.combinevars = [];
    end
    if ~isempty(settings.combinevars)
        settings.changeorder = [];
    end
end