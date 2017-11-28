% script to read and reformat BrainWave-matrix-files
% output is xls-file
%
% [Result,cfg] = lab_convert_bw_matrices(matrix_file,cfg)
%
% written by F. Hatz 2012

function [Result,cfg] = lab_convert_bw_matrices(matrix_file,cfg)

if ~exist('cfg','var')
    cfg = [];
end

if ~exist('matrix_file','var') | ~exist(matrix_file,'file')
    [matrix_file,matrix_filepath]=uigetfile('*.*','Select BrainWave matrix file');
elseif ischar (matrix_file)
    [matrix_file,matrix_filepath] = lab_filename(matrix_file);
end

if ~isnumeric(matrix_file)
    cd(matrix_filepath);
    matrix = lab_read_bw_matrices(fullfile(matrix_filepath,matrix_file));
    [~,~,~,matrix_file] = lab_filename(matrix_file);
else
    matrix = matrix_file;
    clearvars matrix_file;
end

if isempty(matrix) | (isnumeric(matrix) & matrix == 0)
    disp('Abort, no valid matrix')
    Result = [];
    return
end

prompt = matrix.name{1,1};
tmp = union(strfind(prompt,'_'),strfind(prompt,' '));
if ~isempty(tmp)
    if tmp(1) > 1
        subjectname = prompt(1:tmp(1)-1);
    else
        subjectname = [];
    end
    for i = 1:length(tmp)
        tmp2 = num2str(i-1);
        numi = ' ';
        for j = 1:length(tmp2)
            numi = [numi '^' tmp2(j)];
        end
        if i ~= length(tmp)
            subjectname = [subjectname numi ' ' prompt(tmp(i)+1:tmp(i+1)-1)];
        elseif tmp(i) == length(prompt)
            subjectname = [subjectname numi];
        else
            subjectname = [subjectname numi ' ' prompt(tmp(i)+1:end)];
        end
    end
else
    subjectname = prompt;
end
if ~isfield(cfg,'subjectname') | cfg.subjectname < 0
    cfg.subjectname = 0;
end
clearvars prompt tmp tmp2 numi
Prompt = {};
Formats = {};
Prompt(end+1,:) = {'Number of underscores in subject name',''};
Formats(end+1,1).type = 'text';
Prompt(end+1,:) = {['(' subjectname ')'],'subjectname'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'integer';
Formats(end,1).size = 40;
Formats(end,1).limits = [-1 20];
[cfg,Cancelled] = inputsdlg(Prompt,'Convert matrixes',Formats,cfg);
if Cancelled == 1
    return
end
clearvars subjectname Formats Prompt Cancelled i j

namelog = [];
subjectnr = 0;
for i = 1:size(matrix.name,2)
    tmp = strfind(matrix.name{1,i},'_');
    if cfg.subjectname < length(tmp)
        matrix.name{1,i} = matrix.name{1,i}(1:tmp(cfg.subjectname+1)-1);
    end
    if isempty(namelog) | ~strcmp(namelog,matrix.name{1,i})
        subjectnr = subjectnr+1;
        subjecttrials{subjectnr} = i;
        subjectname{subjectnr} = matrix.name{1,i};
        namelog = matrix.name{1,i};
    else
        subjecttrials{subjectnr} = [subjecttrials{subjectnr} i];
    end
end
clearvars namelog i tmp

for i = 1:size(subjectname,2)
    tmp = regexp(subjectname{1,i},'\d');
    if length(tmp) == length(subjectname{1,i})
        subjectname{1,i} = ['P_' subjectname{1,i}];
    end
end

if exist('matrix_file','var')
    warning off
    mkdir(fullfile(matrix_filepath,['M_' matrix_file]));
    warning on
    matrix_filepath = fullfile(matrix_filepath,['M_' matrix_file]);
end

Result.matrix = zeros(size(matrix.matrix,1),size(matrix.matrix,2),subjectnr);
for i = 1:subjectnr
    Result.matrix(:,:,i) = mean(matrix.matrix(:,:,subjecttrials{i}),3);
    Result.name{i} = subjectname{i};
    if exist('matrix_filepath','var')
        dlmwrite(fullfile(matrix_filepath,[subjectname{i} '_matrix.txt']), ...
            Result.matrix(:,:,i),'delimiter','\t','precision', 6);
    end
end
if exist('matrix_filepath','var')
    dlmwrite(fullfile(matrix_filepath,'AVG_matrix.txt'), ...
        mean(Result.matrix,3),'delimiter','\t','precision', 6);
end

button = questdlg('Split by Outcomes-file?','Split by Outcomes-file','Cancel','Yes','No','No');
if strcmp(button,'Yes')
    disp ('Read outcomes file')
    [results_file,results_filepath]=uigetfile('*.xls;*.xlsx','Select file with outcomes');
    if ~results_file == 0 & exist(fullfile(results_filepath,results_file),'file')
        if ispc
            [~,~,results] = xlsread(fullfile(results_filepath,results_file));
        else
            [~,~,results] = xlsread(fullfile(results_filepath,results_file),1,'','basic');
        end
    else
        results = [];
    end
    if ~isempty(results)
        if length(intersect(subjectname,results(2:end,1))) > 2
            results = results';
        elseif length(intersect(subjectname,results(1,2:end))) < 2
            disp('Data and Results not matching')
            return
        end
        subjectsresults = results(1,2:end);
        subjects = intersect(subjectsresults,subjectname);
        subjects = subjects(:)';
        if length(subjects) == length(subjectsresults)
            subjects = subjectsresults;
        end
        for i = 1:size(subjects,2)
            tmp = find(~cellfun(@isempty,strfind(subjectname,subjects{1,i})));
            subjectsnr(1,i) = tmp(1,1);
            tmp = find(~cellfun(@isempty,strfind(subjectsresults,subjects{1,i})));
            subjectsnr(2,i) = tmp(1,1);
        end
        clearvars tmp
        disp ('Ask for outcomes')
        strlist = results(2:end,1);
        outcomes = listdlg('PromptString','Select outcomes:','SelectionMode', ...
            'single','ListString',strlist);
        clearvars strlist
        if ~isempty(outcomes)
            resultstmp = cell2mat(results(outcomes+1,subjectsnr(2,:)+1));
            tmp = unique(resultstmp);
            for i = 1:length(tmp)
                Resultout = mean(Result.matrix(:,:,subjectsnr(1,resultstmp == tmp(i))),3);
                dlmwrite(fullfile(matrix_filepath,[num2str(tmp(i)) '_matrix.txt']), ...
                    Resultout,'delimiter','\t','precision', 6);
                eval(['Result.matrix' num2str(tmp(i)) ' = Resultout;']);
                cfg.Output_filepath = matrix_filepath;
                cfg.Output_file = [matrix_file '_' num2str(tmp(i))];
                if ~exist('settings','var')
                    settings = lab_set_graphanalysis;
                end
                if isfield(settings,'GRAPH') & ~isempty(settings.GRAPH)
                    [Resultgraph,settings,cfg] = lab_graphanalysis(Resultout,[],settings,cfg);
                    eval(['Result.graph' num2str(tmp(i)) ' = Resultgraph;']);
                end
                clearvars Resultout Resultgraph
            end
        end    
    end
end


