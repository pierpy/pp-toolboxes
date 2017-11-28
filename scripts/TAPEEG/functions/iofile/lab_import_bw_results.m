% Script to import BrainWave results, output is a cell-array ('xlsout') or
% a structure ('Result'). Dialogs will be shown to ask, if trials of single
% subjects (identification by name) should be averaged and to select
% variables needed. For cell-array you have to decide between single and
% mean values
%
% written by F. Hatz 2012

function [xlsout,Result,cfg] = lab_import_bw_results(filename,cfg,xlsmode,selection)

xlsout = [];
Result = [];
showdiag = 0;
if ~exist('cfg','var')
    cfg = [];
end

if ~exist('filename','var') | isempty(filename)
    [filename,filepath]=uigetfile('*.*','Select BrainWave results-file');
    if filename == 0
        return
    end
    cd(filepath);
else
    tmp=strfind(filename,filesep);
    if ~isempty(tmp)
        filepath=filename(1:tmp(end));
        filename=filename(tmp(end)+1:end);
    else
        filepath=[];
    end
    clearvars tmp;
end


Result = lab_read_bw_results(fullfile(filepath,filename));
if isempty(Result)
    return
end

% Ask for settings subject info
if ~exist('cfg','var') | ~isfield(cfg,'subjectname')
    prompt = Result.subjects{1};
    tmp = union(strfind(prompt,'_'),strfind(prompt,' '));
    if size(tmp,1) > 1
        tmp = tmp';
    end
    if ~isempty(tmp)
        if tmp(1) > 1
            subjectname = prompt(1:tmp(1)-1);
        else
            subjectname = [];
        end
        for i = 1:size(tmp,2)
            tmp2 = num2str(i-1);
            numi = ' ';
            for j = 1:length(tmp2)
                numi = [numi '^' tmp2(j)]; %#ok<AGROW>
            end
            if i ~= size(tmp,2)
                subjectname = [subjectname numi ' ' prompt(tmp(i)+1:tmp(i+1)-1)]; %#ok<AGROW>
            elseif tmp(i) == length(prompt)
                subjectname = [subjectname numi]; %#ok<AGROW>
            else
                subjectname = [subjectname numi ' ' prompt(tmp(i)+1:end)]; %#ok<AGROW>
            end
        end
    else
        subjectname = prompt;
    end
    cfg.subjectname = 0;
    clearvars prompt tmp tmp2 numi
    showdiag = 1;
end
if ~isfield(cfg,'AVGsubject')
    cfg.AVGsubject = true;
    showdiag = 1;
end
if ~isfield(cfg,'meansingle')
    if exist('xlsmode','var') & ischar(xlsmode)
        if strcmp(xlsmode,'Single') | strcmp(xlsmode,'Mean')
            cfg.meansingle = xlsmode;
        else
            cfg.meansingle = 'No';
        end
    else
        cfg.meansingle = 'Single';
    end
    showdiag = 1;
end
if ~isfield(cfg,'AVGpowerabs')
    if isfield(Result,'readingerror') & Result.readingerror == true
        cfg.AVGpowerabs = false;
    else
        cfg.AVGpowerabs = true;
    end
end

if showdiag == 1
    Prompt = {};
    Formats = {};
    Prompt(end+1,:) = {'Number of underscores in subject name',''};
    Formats(end+1,1).type = 'text';
    
    Prompt(end+1,:) = {['(' subjectname ')'],'subjectname'};
    Formats(end+1,1).type = 'edit';
    Formats(end,1).format = 'integer';
    Formats(end,1).size = 40;
    Formats(end,1).limits = [-1 20];
    
    Prompt(end+1,:) = {'',''};
    Formats(end+1,1).type = 'text';
    
    Prompt(end+1,:) = {'Average trials per subject','AVGsubject'};
    Formats(end+1,1).type = 'check';
    Formats(end,1).size = [-1 -1];
    
    Prompt(end+1,:) = {'Average power with absolute values','AVGpowerabs'};
    Formats(end+1,1).type = 'check';
    Formats(end,1).size = [-1 -1];
    
    Prompt(end+1,:) = {'',''};
    Formats(end+1,1).type = 'text';
    
    Prompt(end+1,:) = {'Mean or single values','meansingle'};
    Formats(end+1,1).type = 'list';
    Formats(end,1).style = 'popupmenu';
    Formats(end,1).format = 'input';
    Formats(end,1).items = {'Single','Mean','No'};
    
    [cfg,Cancelled] = inputsdlg(Prompt,'Import BrainWave',Formats,cfg);
    if Cancelled == 1
        return
    else
        pause(0.2);
    end
end

subjectsnames = {};
for j = 1:size(Result.subjects,2)
    tmp = strfind(Result.subjects{1,j},'_');
    if size(tmp,2) > cfg.subjectname
        subjectstmp = Result.subjects{1,j}(1:tmp(cfg.subjectname + 1)-1);
    elseif ~isempty(tmp)
        subjectstmp = Result.subjects{1,j}(1:tmp(end)-1);
    else
        subjectstmp = Result.subjects{1,j};
    end
    Nsubject = find(strcmp(subjectsnames,subjectstmp));
    if isempty(Nsubject)
        subjectsnames{1,end+1} = subjectstmp; %#ok<AGROW>
        Nsubject = size(subjectsnames,2);
        subjectstrials{Nsubject} = []; %#ok<AGROW>
    end
    subjectstrials{Nsubject} = [subjectstrials{Nsubject} j]; %#ok<AGROW>
    clearvars tmp subjectstmp
end

for i = 1:size(Result.data,3)
    if find(isnan(Result.data(:,:,i)),1)
        resultsvalid(1,i) = 0; %#ok<AGROW>
    else
        resultsvalid(1,i) = 1; %#ok<AGROW>
    end
end
resultsvalid = find(resultsvalid);
for i = 1:size(Result.datamean,1)
    if find(isnan(Result.datamean(i,:)),2)
        resultsvalidmean(1,i) = 0; %#ok<AGROW>
    else
        resultsvalidmean(1,i) = 1; %#ok<AGROW>
    end
end
resultsvalidmean = find(resultsvalidmean);

if size(subjectstrials,2) < size(Result.subjects,2)
    if cfg.AVGsubject == true
        if strcmp(Result.results{1,1},'total_power') & cfg.AVGpowerabs == true
            disp('    calculate absolute power for average')
            Result.data(:,:,2:7) = Result.data(:,:,2:7) .* repmat(Result.data(:,:,1),[1 1 6]);
        end
        for i = 1:size(subjectstrials,2)
            datatmp(:,i,:) = nanmean(Result.data(:,subjectstrials{i},:),2); %#ok<AGROW>
            datameantmp(:,i) = nanmean(Result.datamean(:,subjectstrials{i}),2); %#ok<AGROW>
        end
        if strcmp(Result.results{1,1},'total_power') & cfg.AVGpowerabs == true
            disp('    calculate relative power after average')
            datatmp(:,:,1) = sum(datatmp(:,:,2:7),3);
            datatmp(:,:,2:7) = datatmp(:,:,2:7) ./ repmat(datatmp(:,:,1),[1 1 6]);
            datameantmp(1:7,:) = mean(permute(datatmp(:,:,1:7),[3 2 1]),3);
        end
        Result.data = datatmp;
        Result.datamean = datameantmp;
        Result.subjects = subjectsnames;
        clearvars datatmp datameantmp i
    end
end

for i = 1:size(Result.subjects,2)
    tmp = regexp(Result.subjects{1,i},'\d');
    if length(tmp) == length(Result.subjects{1,i})
        Result.subjects{1,i} = ['P_' Result.subjects{1,i}];
    end
end

[~,~,~,filenameS] = lab_filename(filename);
% Write mean data to xls
if ~isempty(Result.resultsmean)
    for i = 1:size(Result.resultsmean,1)
        Result.resultsmean{i,1} = regexprep(Result.resultsmean{i,1},'_','');
    end
    xlsoutM = {''};
    xlsoutM(2:size(Result.resultsmean,1)+1,1) = Result.resultsmean;
    xlsoutM(1,2:size(Result.subjects,2)+1) = Result.subjects;
    xlsoutM(2:end,2:end) = num2cell(Result.datamean);
    lab_write_xls(fullfile(filepath,[filenameS '_mean.xlsx']),xlsoutM);
end

% Write single data to xls
if ~isempty(Result.results)
    numbervars = size(Result.data,1);
    for i = 1:size(Result.results,1)
        Result.results{i,1} = regexprep(Result.results{i,1},'_','');
    end
    xlsoutS = {''};
    resultsall = [];
    dataall = [];
    for i = 1:size(Result.results,1);
        for j = 1:numbervars
            resultsalltmp{j,1} = [Result.results{i} '_ch' num2str(j)]; %#ok<AGROW>
        end
        dataall = cat(1,dataall,Result.data(:,:,i));
        resultsall = cat(1,resultsall,resultsalltmp);
        clearvars resultsalltmp
    end
    xlsoutS = [xlsoutS Result.subjects];
    xlsoutS(2:size(resultsall,1)+1,1) = resultsall;
    xlsoutS(2:end,2:end) = num2cell(dataall);
    lab_write_xls(fullfile(filepath,[filenameS '.xlsx']),xlsoutS);
end

% select results for full-data
if ~exist('xlsoutM','var') & ~exist('xlsoutS','var')
    cfg.meansingle = 'No';
    xlsout = [];
elseif ~exist('xlsoutM','var') & strcmp(cfg.meansingle,'Mean')
    cfg.meansingle = 'Single';
elseif ~exist('xlsoutS','var') & strcmp(cfg.meansingle,'Single')
    cfg.meansingle = 'Mean';
end
if strcmp(cfg.meansingle,'No')
    xlsout = [];
elseif strcmp(cfg.meansingle,'Single')
    xlsout = xlsoutS(1,:);
    if exist('selection','var')
        if ischar(selection) & strcmp(selection,'all')
            cfg.selection = resultsvalid;
        elseif isnumeric(selection) & max(selection) <= max(resultsvalid)
            cfg.selection = selection;
        end
    end
    if ~exist('cfg','var') | ~isfield(cfg,'selection')
        strdefault = resultsvalid;
        cfg.selection = listdlg('PromptString','Select results:','SelectionMode', ...
            'multiple','ListString',Result.results,'InitialValue',strdefault);
        clearvars strdefault
    end
    if ~isempty(cfg.selection)
        for i = 1:size(cfg.selection,2)
            xlsout = cat(1,xlsout,xlsoutS((cfg.selection(i)-1)*numbervars+2:cfg.selection(i)*numbervars +1,:));
        end
        cfg.clustervars = numbervars;
        xlsout{1,1} = ['C' num2str(numbervars) ' R0'];
    else
        xlsout = [];
    end
elseif strcmp(cfg.meansingle,'Mean')
    xlsout = xlsoutS(1,:);
    if exist('selection','var')
        if ischar(selection) & strcmp(selection,'all')
            cfg.selection = resultsvalidmean;
        elseif isnumeric(selection) & max(selection) <= max(resultsvalidmean)
            cfg.selection = selection;
        end
    end
    if ~exist('cfg','var') | ~isfield(cfg,'selection')
        strdefault = resultsvalidmean;
        cfg.selection = listdlg('PromptString','Select results:','SelectionMode', ...
            'multiple','ListString',Result.resultsmean,'InitialValue',strdefault);
        clearvars strdefault
    end
    if ~isempty(cfg.selection)
        xlsout = cat(1,xlsout,xlsoutM(cfg.selection+1,:));
        cfg.clustervars = 1;
        xlsout{1,1} = 'C1 R0';
    else
        xlsout = [];
    end
end

return