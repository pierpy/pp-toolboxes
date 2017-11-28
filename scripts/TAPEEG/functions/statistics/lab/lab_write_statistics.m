% Write results of statistics
%
% filename = lab_write_statistics(Rstat,header,cfg)

function filename = lab_write_statistics(Rstat,header,cfg)

nameresults = '';
for i = 1:length(header.result)
    header.result{1,i} = regexprep(header.result{1,i},{':',';','/','\'},'');
    nameresults = [nameresults '_' header.result{1,i}]; %#ok<AGROW>
end
if length(nameresults) > 20
    nameresults = nameresults(1,1:20);
end

filename = fullfile(header.path,[header.file nameresults '.xlsx']);
if ~exist(header.path,'dir')
    mkdir(header.path);
end
if exist(filename,'file')
    delete(filename);
end

%initialize verbose file
fid=fopen([filename(1:end-4) '.vrb'],'w');

for Nresult = 1:cfg.resultvars
    if isfield(cfg,'write_meanstd') & length(cfg.write_meanstd) >= Nresult
        write_meanstd = cfg.write_meanstd(Nresult);
    else
        write_meanstd = false;
    end
    if isfield(cfg,'clustervars2') & cfg.clustervars2 > 1
        Result = {['C' num2str(cfg.clustervars) ' R0 V' num2str(cfg.clustervars2)]};
    elseif isfield(cfg,'clustervars')
        Result = {['C' num2str(cfg.clustervars) ' R0']};
    else
        Result = {'Results'};
    end
    Result = cat(2,Result,header.vars(:)');
    if isfield(Rstat{Nresult},'T')
        Result(end+1,1) = {'T'}; %#ok<AGROW>
        Result(end,2:end) = num2cell(Rstat{Nresult}.T);
    end
    if isfield(Rstat{Nresult},'F')
        Result(end+1,1) = {'F'}'; %#ok<AGROW>
        Result(end,2:end) = num2cell(Rstat{Nresult}.F);
    end
    if isfield(Rstat{Nresult},'R')
        Result(end+1,1) = {'R'}; %#ok<AGROW>
        Result(end,2:end) = num2cell(Rstat{Nresult}.R);
    end
    if isfield(Rstat{Nresult},'Z')
        Result(end+1,1) = {'Z'}; %#ok<AGROW>
        Result(end,2:end) = num2cell(Rstat{Nresult}.Z);
    end
    Result(end+1,1) = {'p'}; %#ok<AGROW>
    Result(end,2:end) = num2cell(Rstat{Nresult}.p);
    if isfield(Rstat{Nresult},'p_fdr')
        Result(end+1,1) = {'p_FDR'}; %#ok<AGROW>
        Result(end,2:end) = num2cell(Rstat{Nresult}.p_fdr);
    end
    if isfield(Rstat{Nresult},'mxp')
        Result(end+1,1) = {'mxp'}; %#ok<AGROW>
        Result(end,2:end) = num2cell(Rstat{Nresult}.mxp);
    end
    if isfield(Rstat{Nresult},'mxpV')
        Result(end+1,1) = {'mxpV'}; %#ok<AGROW>
        Result(end,2:end) = num2cell(Rstat{Nresult}.mxpV);
    end
    if isfield(Rstat{Nresult},'mxpC')
        Result(end+1,1) = {'mxpC'}; %#ok<AGROW>
        Result(end,2:end) = num2cell(Rstat{Nresult}.mxpC);
    end
    if isfield(Rstat{Nresult},'mxpS')
        Result(end+1,1) = {'mxpS'}; %#ok<AGROW>
        Result(end,2:end) = num2cell(Rstat{Nresult}.mxpS);
    end
    if isfield(Rstat{Nresult},'df')
        Result(end+1,1) = {'df'}'; %#ok<AGROW>
        Result(end,2:end) = num2cell(Rstat{Nresult}.df);
    end
    if isfield(Rstat{Nresult},'mean') & write_meanstd == true
        for j = 1:length(Rstat{Nresult}.mean)
            Result{end+1,1} = ['Mean_' Rstat{Nresult}.group{j}]; %#ok<AGROW>
            Result(end,2:end) = num2cell(Rstat{Nresult}.mean{j});
            Result{end+1,1} = ['Std_' Rstat{Nresult}.group{j}]; %#ok<AGROW>
            Result(end,2:end) = num2cell(Rstat{Nresult}.std{j});
            Result{end+1,1} = ['Median_' Rstat{Nresult}.group{j}]; %#ok<AGROW>
            Result(end,2:end) = num2cell(Rstat{Nresult}.median{j});
            Result{end+1,1} = ['LowerQuartile_' Rstat{Nresult}.group{j}]; %#ok<AGROW>
            Result(end,2:end) = num2cell(Rstat{Nresult}.percent25{j});
            Result{end+1,1} = ['UpperQuartile_' Rstat{Nresult}.group{j}]; %#ok<AGROW>
            Result(end,2:end) = num2cell(Rstat{Nresult}.percent75{j});
            Result{end+1,1} = ['Min_' Rstat{Nresult}.group{j}]; %#ok<AGROW>
            Result(end,2:end) = num2cell(Rstat{Nresult}.min{j});
            Result{end+1,1} = ['Max_' Rstat{Nresult}.group{j}]; %#ok<AGROW>
            Result(end,2:end) = num2cell(Rstat{Nresult}.max{j});
        end
    end
    warning off %#ok<WNOFF>
    lab_write_xls(filename, Result',header.result{1,Nresult});
    if isfield(Rstat{Nresult},'stats') & ~isempty(Rstat{Nresult}.stats)
        ResultLR = [{'LogRegr'} header.vars];
        ResultLR = cat(1,ResultLR,Rstat{Nresult}.stats);
        VarsNr = (size(ResultLR,1) - 3) / 3;
        ResultLR{4,1} = 'B0';
        ResultLR{4+VarsNr,1} = 'T0';
        ResultLR{4+VarsNr*2,1} = 'pvalue0';
        for i = 1:VarsNr-2
            ResultLR{(4+i),1} = ['B_' header.factors{1,i}];
            ResultLR{(4+i+VarsNr),1} = ['T_' header.factors{1,i}];
            ResultLR{(4+i+VarsNr*2),1} = ['pvalue_' header.factors{1,i}];
        end
        ResultLR{end-2*VarsNr,1} = ['B_' header.result{1,Nresult}];
        ResultLR{end-VarsNr,1} = ['T_' header.result{1,Nresult}];
        ResultLR{end,1} = ['pvalue_' header.result{1,Nresult}];
        lab_write_xls(filename, ResultLR',['LogRegr_' header.result{1,Nresult}]);
        clearvars ResultLR
    end
    warning on %#ok<WNON>
    
    % write verbose info
    fprintf(fid,['OUTCOME: ' header.result{1,Nresult} '\n']);
    if cfg.permutations(Nresult) > 1
        fprintf(fid,'Permutation\n');
        fprintf(fid,datestr(now,0));
        fprintf(fid,'\n\n');
        fprintf(fid,'Number of permuations\n');
        fprintf(fid,num2str(cfg.permutations(Nresult)));
        fprintf(fid,'\n\n');
    end
    if cfg.clustervars > 1
        fprintf(fid,'number of Variables\n');
        fprintf(fid,num2str(cfg.clustervars));
        if isfield(cfg,'clustervarsI')
            fprintf(fid,[' (' num2str(cfg.clustervarsI) ' included)']);
        end
        fprintf(fid,'\n\n');
        fprintf(fid,'Variables\n');
        fprintf(fid,sprintf('%s ',header.vars{:}));
        fprintf(fid,'\n\n');
    end
    fprintf(fid,'Number of Measures\n');
    fprintf(fid,num2str(cfg.numclusters));
    fprintf(fid,'\n\n');
    if isfield(header,'measures')
        fprintf(fid,'Measures\n');
        fprintf(fid,sprintf('%s ',header.measures{:}));
        fprintf(fid,'\n\n');
    end
    fprintf(fid,'Method\n');
    fprintf(fid,cfg.method{Nresult,1});
    fprintf(fid,'\n\n');
    fprintf(fid,'Outcome\n');
    fprintf(fid,header.result{1,Nresult});
    fprintf(fid,'\n\n');
    if isfield(cfg,'factorsvars')
        fprintf(fid,'Number of factors for generalized linear model regression\n');
        fprintf(fid,num2str(cfg.factorvars));
        if ~isempty(header.factors)
            fprintf(fid,sprintf('%s ',header.factors{:}));
        end
        fprintf(fid,'\n\n');
    end
    fprintf(fid,'Subjects\n');
    fprintf(fid,sprintf('%s ',header.subjects{:}));
    if isfield(cfg,'filelist') & ~isempty(cfg.filelist)
        fprintf(fid,'\n\n');
        fprintf(fid,'Input files\n');
        fprintf(fid,sprintf('%s ',cfg.filelist{:}));
    end
    fprintf(fid,'\n\n');
end

% close verbose file
fclose(fid);

% save matlab container file with results
save([filename(1:end-4) '.mat'],'Rstat','header','cfg');
