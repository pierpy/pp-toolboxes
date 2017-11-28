% Repeated measures Anova
%
% [data,header,result,factors,cfg] = lab_read_statistics(cfg,doresults,dofactors,doclusters,dowriting,doselection)
%
% Read xls-file with results for statistics
%
%    doresults       1 = select outcomes for statistics
%    dorfactors      1 = select factors for permutation-statistics
%    doclusters      1 = for every measure multiple variables are stored in xls-file
%                              measure1_variable1
%                              measure1_variable2
%                              measure1_variable3
%                              measure2_variable1
%                              measure2_variable2
%                              ...
%    dowriting      1 = selected data for statistics is stored in a new xls-file
%    doselection    1 = selected measures,outcomes,factors and preprocessings
%
% written by F. Hatz 2013

function [data,header,result,factors,cfg] = lab_read_statistics(cfg,doresults,dofactors,doclusters,dowriting,doselection)

data = [];
result = [];
factors = [];
header = [];

if ~exist('doselection','var') | isempty(doselection)
    doselection = 1;
end
if ~exist('dowriting','var') | isempty(dowriting)
    dowriting = 1;
end
if ~exist('doclusters','var') | isempty(doclusters)
    doclusters = 1;
end
if ~exist('dofactors','var') | isempty(dofactors)
    dofactors = 1;
end
if ~exist('doresults','var') | isempty(doresults)
    doresults = 1;
end
if ~exist('cfg','var')
    cfg = [];
end


if ~isempty(cfg) & iscell(cfg)
    datainput = cfg;
    cfg = [];
else
    if ischar(cfg) & exist(cfg,'file')
        [data_file,data_filepath] = lab_filename(cfg);
        cfg = [];
        readdata = 0;
        datainput = [];
    else
        if isfield(cfg,'selectmode') & (strcmp(cfg.selectmode,'Multiple') | strcmp(cfg.selectmode,'Single'))
            selectmode = cfg.selectmode;
            cfg = rmfield(cfg,'selectmode');
        else
            disp('Multiple files or single file')
            selectmode = questdlg('Multiple or single input file(s)?','Single/Multiple','Cancel','Multiple','Single','Single');
        end
        if strcmp(selectmode,'Multiple')
            disp ('Read file')
            [data_file,data_filepath]=uigetfile('*.*','Select file Nr: 1');
            readdata = 1;
        elseif strcmp(selectmode,'Single')
            disp ('Read file')
            [data_file,data_filepath]=uigetfile('*.*','Select File');
            readdata = 0;
        else
            return
        end
        datainput = [];
    end
    if isempty(data_file) | data_file == 0
        return
    end
end

if isempty(datainput)
    cd(data_filepath);
    if strcmp(data_file(end-2:end),'xls') | strcmp(data_file(end-3:end),'xlsx') | strcmp(data_file(end-2:end),'csv')
        if strcmp(data_file(end-2:end),'csv')
            datainput = lab_read_csv(fullfile(data_filepath,data_file));
        else
            datainput = lab_read_xls(fullfile(data_filepath,data_file));
        end
        disp('Subjects in row or column')
        Mtranspose = questdlg('Are subjects names in first row or column?','Subjects column/row','Cancel','Column','Row','Row');
        if isempty(Mtranspose) | strcmp(Mtranspose,'Cancel')
            return
        elseif strcmp(Mtranspose,'Column')
            datainput = datainput';
            datainput{1,1} = '';
        end
        datainput = lab_correctheader(datainput);
        [datainput,cfg] = lab_getstructure(datainput,cfg);
        datainput = lab_set_NaN(datainput);
        if doclusters == 0
            cfg.clustervars = 1;
            if isfield(cfg,'clustervars2')
                cfg = rmfield(cfg,'clustervars2');
            end
        end
        clustervars = cfg.clustervars;
        numresults = cfg.numresults;
        cfg.filelist = cellstr(data_file);
        grouptmp = cellstr('group');
        grouptmp(1,2:size(datainput,2)) = repmat(num2cell(1),1,size(datainput,2)-1);
        postfix.file1 = lab_correct_freq(data_file);
        if isfield(cfg,'collectfiles') & cfg.collectfiles == 1
            postfix.mode = 'measures (post)';
            cfg = rmfield(cfg,'collectfiles');
        end
        cfg.numfiles = 1;
        while readdata > 0
            [data2_file,data2_filepath]=uigetfile('*.xls;*.xlsx',['Select file nr: ' num2str(readdata+1)]);
            if data2_file ~= 0
                if strcmp(data2_file(end-2:end),'csv')
                    datainput2 = lab_read_csv(fullfile(data2_filepath,data2_file));
                else
                    datainput2 = lab_read_xls(fullfile(data2_filepath,data2_file));
                end
                if strcmp(Mtranspose,'Column')
                    datainput2 = datainput2';
                    datainput2{1,1} = '';
                end
                datainput2 = lab_correctheader(datainput2);
                [datainput2,cfg2] = lab_getstructure(datainput2);
                datainput2 = lab_set_NaN(datainput2);
                if doclusters == 0
                    cfg2.clustervars = 1;
                    if isfield(cfg2,'clustervars2')
                        cfg2 = rmfield(cfg2,'clustervars2');
                    end
                end
                if ~isempty(cfg2.clustervars) & ~isempty(clustervars) & cfg2.clustervars ~= clustervars
                    disp('Exclude last file, number of variables per measure not matching')
                    readdata = 0;
                    datainput2 = [];
                end
                if ~isempty(cfg2.numresults) & ~isempty(numresults) & cfg2.numresults ~= numresults
                    disp('Exclude last file, number of results not matching')
                    readdata = 0;
                    datainput2 = [];
                end
                postfix.file2 = lab_correct_freq(data2_file);
                [datainput,grouptmp,cfg,postfix,skipprocessing] = lab_read_statistics_combine(datainput,datainput2,grouptmp,cfg,postfix);
                if skipprocessing == 0
                    readdata = readdata + 1;
                    group = grouptmp;
                    cfg.filelist = [cfg.filelist cellstr(data2_file)];
                    cd(data2_filepath);
                else
                    cfg.numfiles = readdata;
                    readdata = 0;
                end
            else
                cfg.numfiles = readdata;
                readdata = 0;
            end
            if exist('data2_filepath','var') & ischar(data2_filepath)
                cd(data2_filepath);
            end
            clearvars data2_file data2_filepath
        end
        clearvars grouptmp
    else
        [datainput,~,cfg] = lab_import_bw_results(fullfile(data_filepath,data_file),cfg);
        cfg.filelist = cellstr(data_file);
        grouptmp{1,1} = 'group';
        grouptmp(1,2:size(datainput,2)) = repmat(num2cell(1),1,size(datainput,2)-1);
        postfix.file1 = lab_correct_freq(data_file);
        cfg.numfiles = 1;
        while readdata > 0
            [data2_file,data2_filepath]=uigetfile('*.*',['Select file nr: ' num2str(readdata+1)]);
            if data2_file ~= 0
                try
                    [datainput2,~,cfg] = lab_import_bw_results(fullfile(data2_filepath,data2_file),cfg);
                catch %#ok<CTCH>
                    datainput2 = [];
                end
                postfix.file2 = lab_correct_freq(data2_file);
                [datainput,grouptmp,cfg,postfix,skipprocessing] = lab_read_statistics_combine(datainput,datainput2,grouptmp,cfg,postfix);
                if skipprocessing == 0
                    readdata = readdata + 1;
                    group = grouptmp;
                    cfg.filelist = [cfg.filelist cellstr(data2_file)];
                else
                    cfg.numfiles = readdata;
                    readdata = 0;
                end
            else
                cfg.numfiles = readdata;
                readdata = 0;
            end
            if exist('data2_filepath','var') & ischar(data2_filepath)
                cd(data2_filepath);
            end
            clearvars data2_file data2_filepath
        end
        clearvars grouptmp
        dowriting = 1;
    end
else
    disp('Subjects in row or column')
    Mtranspose = questdlg('Are subjects names in first row or column?','Subjects column/row','Cancel','Column','Row','Row');
    if isempty(Mtranspose) | strcmp(Mtranspose,'Cancel')
        return
    elseif strcmp(Mtranspose,'Column')
        datainput = datainput';
        datainput{1,1} = '';
    end
    datainput = lab_correctheader(datainput);
    [datainput,cfg] = lab_getstructure(datainput,cfg);
    if doclusters == 0
        cfg.clustervars = 1;
        if isfield(cfg,'clustervars2')
            cfg = rmfield(cfg,'clustervars2');
        end
    end
    cfg.filelist = cellstr('');
    grouptmp = cellstr('group');
    grouptmp(1,2:size(datainput,2)) = repmat(num2cell(1),1,size(datainput,2)-1);
    postfix.file1 = '';
    cfg.numfiles = 1;
end

if doclusters == 0
    cfg.clustervars = 1;
    if isfield(cfg,'clustervars2')
        cfg = rmfield(cfg,'clustervars2');
    end
end

if doresults > 0
    strlist = {'Included','Select file','Select manually'};
    if exist('group','var')
        strlist = [strlist {'Groups by files'}];
    end
    if ~isempty(cfg.clustervars) & cfg.clustervars > 1
        strlist = [strlist {'Groups by clusters'}];
    end
    strlist = [strlist {'None'}];
    outcomesmode = listdlg('PromptString','Select outcomes:','SelectionMode','single','ListString',strlist);
    pause(0.2);
    if isempty(outcomesmode)
        doresults = 0;
        outcomesmode = 'None';
    else
        outcomesmode = strlist(outcomesmode);
        if strcmp(outcomesmode,'None');
            doresults = 0;
        end
    end
else
    doresults = 0;
    outcomesmode = 'None';
end
skiprepeated = false;
if strcmp(outcomesmode,'Included');
    if isempty(cfg.clustervars)
        cfg.clustervars = 1;
    end
    if isempty(cfg.numresults) & cfg.clustervars > 1
        cfg.numresults = mod(size(datainput,1)-1,cfg.clustervars);
    end
elseif strcmp(outcomesmode,'Select file');
    [datainput,cfg] = lab_match_subjects(datainput,cfg);
elseif strcmp(outcomesmode,'Select manually');
    strlist = datainput(1,2:end);
    numbertrack = 1:size(strlist,2);
    group = 1;
    while ~isempty(strlist)
        groupdef = listdlg('PromptString',['Select group ' num2str(group)],'SelectionMode','multiple','ListString',strlist);
        if ~isempty(groupdef)
            tmp = setdiff(1:size(strlist,2),groupdef);
            strlist = strlist(1,tmp);
            groupdef = numbertrack(1,groupdef);
            grouptmp(1,groupdef) = group;
            numbertrack = numbertrack(1,tmp);
            clearvars tmp
            group = group + 1;
        else
            strlist = [];
        end
    end
    group = [{'group'} num2cell(grouptmp)];
    datainput = cat(1,datainput,group);
    if isfield(cfg,'numresults')
        cfg.numresults = cfg.numresults + 1;
    else
        cfg.numresults = 1;
    end
elseif strcmp(outcomesmode,'Groups by files');
    group{1,1} = 'FileNr';
    datainput = cat(1,datainput,group);
    if isfield(cfg,'numresults')
        cfg.numresults = cfg.numresults + 1;
    else
        cfg.numresults = 1;
    end
elseif strcmp(outcomesmode,'Groups by clusters');
    if ~isfield(cfg,'numresults') | isempty(cfg.numresults)
        cfg.numresults = 0;
    end
    Subjects = datainput(1,2:end);
    Vars = datainput(2:cfg.clustervars+1,1);
    for i = 1:length(Vars)
        tmp = strfind(Vars{i},'_');
        if ~isempty(tmp)
            Vars{i} = Vars{i}(tmp(end)+1:end);
        end
    end
    Measures = datainput(2:cfg.clustervars:end-cfg.numresults,1);
    for i = 1:length(Measures)
        tmp = strfind(Measures{i},'_');
        if ~isempty(tmp)
            Measures{i} = Measures{i}(1:tmp(end)-1);
        end
    end
    clearvars tmp
    if cfg.numresults > 0
        Result = datainput(end-cfg.numresults+1:end,2:end);
        ResultV = datainput(end-cfg.numresults+1:end,1);
    else
        Result = {};
        ResultV = {};
    end
    Data = datainput(2:end-cfg.numresults,2:end);
    Data = reshape(Data,cfg.clustervars,length(Measures),size(Data,2));
    Data = permute(Data,[1 3 2]);
    Data = reshape(Data,size(Data,1),size(Data,2)*size(Data,3));
    Result = repmat(Result,[1 length(Measures)]);
    datainput = cat(1,Data,Result);
    datainput = cat(2,cat(1,Vars,ResultV),datainput);
    headertmp = {' '};
    for i = 1:length(Measures)
        for j = 1:length(Subjects)
            headertmp{1,end+1} = [Measures{i} '_' Subjects{j}]; %#ok<AGROW>
        end
    end
    datainput = cat(1,headertmp,datainput);
    Group = {'Cluster'};
    for i = 1:length(Measures)
        Group = [Group repmat(Measures(i),1,length(Subjects))]; %#ok<AGROW>
    end
    datainput = cat(1,datainput,Group);
    cfg.numresults = cfg.numresults+1;
    if isfield(cfg,'clustervars2') & ~isempty(cfg.clustervars2)
        cfg.clustervars2 = 1;
    end
    cfg.clustervars = 1;
    datainput{1,1} = {['C' num2str(cfg.clustervars) ' R' num2str(cfg.numresults+1)]};
    clearvars Measures Vars Data Result ResultV Subjects i j Group
    skiprepeated = true;
elseif ~strcmp(outcomesmode,'None');
    return
end

if exist('postfix','var') & isfield(postfix,'enable')
    [datainput,cfg] = clustervarsnew(datainput,cfg);
elseif skiprepeated == false
    [datainput,cfg] = repeatedtrials(datainput,cfg);
end

if doresults == -1
    cfg.numresults = 0;
    datainput = datainput(1:end-cfg.numresults,:);
end

if exist('data_file','var')
    [~,~,~,cfg.file] = lab_filename(data_file);
    cfg.path = data_filepath;
    cd(cfg.path);
    if exist(fullfile(cfg.path,[cfg.file '.xls']),'file') | exist(fullfile(cfg.path,[cfg.file '.xlsx']),'file')
        [~,~,~,cfg.file] = lab_filename(cfg.file);
        cfg.file = [cfg.file '_statistics'];
    end
else
    cfg.path = pwd;
    cfg.file = '';
end

if doselection >= 1
    tmp = cellfun(@isempty,datainput);
    datainput(tmp) = {NaN};
    tmp = cellfun(@ischar,datainput);
    tmp(:,1) = false;
    tmp(1,:) = false;
    if isfield(cfg,'numresults') & ~isempty(cfg.numresults) & cfg.numresults > 0
        tmp(end-cfg.numresults+1:end,:) = false;
    end
    datainput(tmp) = {NaN};
    tmp = sum(tmp,2) < (size(tmp,2)-1);
    datainput = datainput(tmp,:);
    datainput{1,1} = ['C' num2str(cfg.clustervars) ' R' num2str(cfg.numresults)];
    [data,header,result,factors,cfg] = lab_preprocess_statistics(datainput,cfg,doresults,dofactors,doselection,dowriting);
    if isempty(data)
        return
    end
else
    xlsout = datainput';
    if isfield(cfg,'clustervars2') & cfg.clustervars2 > 1
        xlsout{1,1} = ['C' num2str(cfg.clustervars) ' R' num2str(cfg.numresults) ' V' num2str(cfg.clustervars2)];
    else
        xlsout{1,1} = ['C' num2str(cfg.clustervars) ' R' num2str(cfg.numresults)];
    end
    data = datainput;
    data{1,1} = xlsout{1,1};
    if isfield(cfg,'filelist') & size(cfg.filelist,2) > 1
        if dowriting == 1 & ~isempty(cfg.file)
            tmp.filename = fullfile(cfg.path,[cfg.file '.xlsx']);
            Prompt(1,1:2) = {'','filename'};
            Formats.type = 'edit';
            Formats.format = 'file';
            Formats.items = {'*.xls;*.xlsx','Excel-file (*.xls/xlsx)'};
            Formats.limits = [1 0]; % use uiputfile
            Formats.size = 300;
            [tmp,Cancelled] = inputsdlg(Prompt,'Save File to',Formats,tmp);
            if Cancelled == 1
                data = [];
                result = [];
                factors = [];
                header = [];
                return
            else
                [~,tmp2,~,tmp1] = lab_filename(tmp.filename);
                cfg.file = tmp1;
                cfg.path = tmp2;
                cd(cfg.path);
                clearvars tmp tmp1 tmp2
            end
        else
            Prompt(1,1:2) = {'New filename','file'};
            Formats.type = 'edit';
            Formats.format = 'text';
            Formats.size = 150;
            [cfg,Cancelled] = inputsdlg(Prompt,'Filename',Formats,cfg);
            if Cancelled == 1
                data = [];
                result = [];
                factors = [];
                header = [];
                return
            end
        end
    end
end
header.file = cfg.file;
header.path = cfg.path;

if dowriting == 1 & ~isempty(cfg.file)
    if ~exist('xlsout','var')
        % store cluster-info and result-info in first field
        if isfield(cfg,'clustervars2') & cfg.clustervars2 > 1
            xlsout{1,1} = ['C' num2str(cfg.clustervars) ' R' num2str(cfg.numresults) ' V' num2str(cfg.clustervars2)];
        else
            xlsout{1,1} = ['C' num2str(cfg.clustervars) ' R' num2str(cfg.numresults)];
        end
        xlsout = [xlsout;header.subjects];
        datatmp = cat(1,header.vars,num2cell(data));
        if ~isempty(factors)
            datatmp = [datatmp cat(1,header.factors,num2cell(factors))];
        end
        if ~isempty(result)
            datatmp = [datatmp cat(1,header.result,num2cell(result))];
        end
        xlsout = [xlsout datatmp];
    end
    if ~isempty(xlsout)
        if ~exist(cfg.path,'dir')
            mkdir(cfg.path);
        end
        lab_write_xls(fullfile(cfg.path,[cfg.file '.xlsx']),xlsout');
    end
end

end



function [datainput,cfg] = clustervarsnew(datainput,cfg)
   if ~isfield(cfg,'numresults') | isempty(cfg.numresults)
       cfg.numresults = 0;
   end
   for i = 2:size(datainput,1)-cfg.numresults
       tmp = strfind(datainput{i,1},'_');
       if ~isempty(tmp)
           tmp2{i-1,1} = datainput{i,1}(tmp(end):end); %#ok<AGROW>
       else
           tmp2{i-1,1} = datainput{i,1}; %#ok<AGROW>
       end
   end
   [~,m,~]=unique(tmp2,'last');
   if max(mod(m,m(1))) == 0
       datatmp = datainput(2:end-cfg.numresults,1:end);
       datatmp = reshape(datatmp,[m(1), size(datatmp,1)/m(1) size(datatmp,2)]);
       datatmp = permute(datatmp,[2 1 3]);
       datatmp = reshape(datatmp,[size(datatmp,1)*size(datatmp,2) size(datatmp,3)]);
       datainput(2:end-cfg.numresults,1:end) = datatmp;
       clearvars datatmp
   end
   for i = 2:size(datainput,1)-cfg.numresults
       tmp = strfind(datainput{i,1},'_');
       if ~isempty(tmp)
           tmp2{i-1,1} = datainput{i,1}(1:tmp(end)-1); %#ok<AGROW>
       else
           tmp2{i-1,1} = datainput{i,1}; %#ok<AGROW>
       end
   end
   [~,m,~]=unique(tmp2,'last');
   m = sort(m);
   n = mod(m,m(1));
   if n == 0 & (length(m) == 1 | (m(1) > 1 & (m(2)-m(1)) == m(1)))
       cfg.clustervars = m(1);
   end
   clearvars tmp tmp2 m n i
end

function [datainput,cfg] = repeatedtrials(datainput,cfg)
   if size(datainput,2) < 3
       return
   end
   for i = 2:size(datainput,2)
       tmp = strfind(datainput{1,i},'_');
       if ~isempty(tmp)
           tmp1{i-1,1} = datainput{1,i}(1:tmp(end)-1); %#ok<AGROW>
           tmp2{i-1,1} = datainput{1,i}(tmp(end)+1:end); %#ok<AGROW>
       else
           tmp1{i-1,1} = datainput{1,i}; %#ok<AGROW>
           tmp2{i-1,1} = datainput{1,i}; %#ok<AGROW>
       end
       clearvars tmp
   end
   if verLessThan('matlab','8')
       [Trials,m,~]=unique(tmp1,'last');
       [Rvariables,n] = unique(tmp2,'last');
   else
       [Trials,m,~]=unique(tmp1,'stable');
       [Rvariables,n] = unique(tmp2,'stable');
   end
   Measures = datainput(2:end,1);
   if length(Trials) <= 1
       return
   end
   if max(mod(m-1,length(Rvariables))) == 0
       mode = 1;
       flag = 1;
   elseif max(mod(n-1,length(Trials))) == 0
       mode = 2;
       flag = 1;
   else
       flag = 0;
   end
   if flag == 0
       return
   end
   if isfield(cfg,'arrangemode')
       arrangemode = cfg.arrangemode;
       cfg = rmfield(cfg,'arrangemode');
   else
       arrangemode = questdlg('Re-arrange repeated subjects/trials?','Re-arrange','No','Yes','No');
   end
   if strcmp(arrangemode,'Yes')
       data = datainput(2:end,2:end);
       if mode == 1
           data = reshape(data,[size(data,1) length(Rvariables) size(data,2)/length(Rvariables)]);
           data = permute(data,[2 1 3]);
       elseif mode == 2
           data = reshape(data,[size(data,1) size(data,2)/length(Rvariables) length(Rvariables)]);
           data = permute(data,[3 1 2]);
       else
           return
       end
       data = reshape(data,size(data,1)*size(data,2),size(data,3));
       MeasuresOut = {};
       cfg.numresults = cfg.numresults * length(Rvariables);
       for i = 1:length(Measures)
           for j = 1:length(Rvariables)
               MeasuresOut{end+1,1} = [Measures{i} '_' Rvariables{j}]; %#ok<AGROW>
           end
       end
       if isfield(cfg,'numresults') & cfg.numresults > 0
           dataout = cellstr(['C' num2str(size(data,3)) '_R' num2str(cfg.numresults)]);
       else
           dataout = cellstr(['C' num2str(size(data,3)) '_R0']);
       end
       dataout = [dataout Trials(:)'];
       dataout = cat(1,dataout,[MeasuresOut data]);
       datainput = dataout;
       if isfield(cfg,'clustervars') & cfg.clustervars > 1
           cfg.clustervars2 = cfg.clustervars;
       end
       cfg.clustervars = length(Rvariables);
   elseif mode == 2
       datatmp = datainput(:,2:end);
       tmp = 1:size(datatmp,2);
       tmp = reshape(tmp,[size(datatmp,2)/length(Rvariables) length(Rvariables)])';
       tmp = tmp(:)';
       datainput(:,2:end) = datatmp(:,tmp);
   end
end