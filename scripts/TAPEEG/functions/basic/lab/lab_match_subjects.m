function [datainput,cfg] = lab_match_subjects(datainput,cfg,results)

disp ('Read outcomes file')
    
if ~exist('results','var')
    [results_file,results_filepath]=uigetfile('*.xls;*.xlsx','Select file with outcomes');
    if ~results_file == 0 & exist(fullfile(results_filepath,results_file),'file')
        results = lab_read_xls(fullfile(results_filepath,results_file));
        results = lab_set_NaN(results);
    else
        results = [];
    end
end
    
if isempty(datainput) | isempty(results) | ~iscell(results)
    return
end
datainput = lab_correctheader(datainput);
results = lab_correctheader(results);
Dsubjects1 = upper(regexprep(datainput(1,2:end),{'-','_',' '},''));
Dsubjects2 = upper(regexprep(datainput(2:end,1),{'-','_',' '},''));
Rsubjects1 = upper(regexprep(results(1,2:end),{'-','_',' '},''));
Rsubjects2 = upper(regexprep(results(2:end,1),{'-','_',' '},''));

if length(intersect(Dsubjects1,Rsubjects1)) > 2
    Dsubjects = Dsubjects1;
    Rsubjects = Rsubjects1;
elseif length(intersect(Dsubjects1,Rsubjects2)) > 2
    Dsubjects = Dsubjects1;
    Rsubjects = Rsubjects2;
    results = results';
elseif length(intersect(Dsubjects2,Rsubjects1)) > 2
    Dsubjects = Dsubjects2;
    Rsubjects = Rsubjects1;
    datainput = datainput';
elseif length(intersect(Dsubjects2,Rsubjects2)) > 2
    Dsubjects = Dsubjects2;
    Rsubjects = Rsubjects2;
    datainput = datainput';
    results = results';
else
    disp('Data and Results not matching')
    return
end
subjects = intersect(Rsubjects,Dsubjects);
subjects = subjects(:)';
if length(subjects) == length(Rsubjects)
    subjects = Rsubjects(:)';
end
for i = 1:length(subjects)
    tmp = find(~cellfun(@isempty,strfind(Dsubjects,subjects{1,i})));
    subjectsnr(1,i) = tmp(1,1); %#ok<AGROW>
    tmp = find(~cellfun(@isempty,strfind(Rsubjects,subjects{1,i})));
    subjectsnr(2,i) = tmp(1,1); %#ok<AGROW>
end
clearvars tmp
tmp = min(cellfun(@isnumeric,results(2:end,2:end)),[],2);
tmp = logical([0;tmp]);
datainput = datainput(:,[1 (subjectsnr(1,:)+1)]);
datainput = cat(1,datainput,results(tmp,[1 subjectsnr(2,:)+1]));
if isfield(cfg,'numresults')
    cfg.numresults = cfg.numresults + sum(tmp);
else
    cfg.numresults = sum(tmp);
end