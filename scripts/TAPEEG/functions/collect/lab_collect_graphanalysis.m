% Collect results from graphanalysis analysis
% Output is xls-files
%
% lab_collect_graphanalysis(searchfolder)
%
% written by F. Hatz 2013

function lab_collect_graphanalysis(cfg)
disp('Collect graphanalysis data of all files')

skipprocessing = 0;
if ~exist('cfg','var')
    cfg = [];
    skipselection = false;
else
    skipselection = true;
end
FilesAll = [];

if ~isfield(cfg,'CollectGraph') | ~isfield(cfg.CollectGraph,'searchfolder')
    [cfg,FilesAll,skipprocessing] = lab_set_collect_graphanalysis(cfg);
    if skipprocessing == 1
        return
    end
end

% turn log-file on
if exist(cfg.CollectGraph.searchfolder,'dir')
    diary(fullfile(cfg.CollectGraph.searchfolder,'CollectGraph.log'));
end

% search files
if isempty(FilesAll)
    FilesAll = lab_collect_graphanalysis_search(cfg,skipselection);
    if skipprocessing == 1
        diary off
        return
    end
end

if isempty(FilesAll)
    disp('no graph results to collect')
    return
end

if isempty(cfg.CollectGraph.outputfolder)
    disp('Abort: specify output folder')
    return
end
warning off %#ok<WNOFF>
mkdir(fullfile(cfg.CollectGraph.searchfolder,cfg.CollectGraph.outputfolder));
Output_path = fullfile(cfg.CollectGraph.searchfolder,cfg.CollectGraph.outputfolder);
warning on %#ok<WNON>

for connectnr = 1:size(FilesAll,2)
    Files = FilesAll(1,connectnr);
    for filenr = 1:length(Files.list)
        skipfile = 0;
        disp(['Collect graphanalysis data: ' lab_filename(Files.list{filenr})])
        MAT = load(Files.list{filenr});
        if ~isfield(MAT,'Result')
            skipfile = 1;
        end
        if skipfile == 0
            if isfield(MAT.Result,'data') & ~isempty(MAT.Result.data)
                numchans = size(MAT.Result.data,1);
            elseif isfield(MAT.Result,'dataAVG') & ~isempty(MAT.Result.dataAVG)
                numchans = size(MAT.Result.dataAVG,1);
            else
                skipfile = 1;
            end
        end
        if skipfile == 0
            if isfield(MAT.Result,'freqband') & ~isempty(MAT.Result.freqband)
                freqband = MAT.Result.freqband;
            else
                tmp = strfind(Files.list{filenr},'_F');
                if ~isempty(tmp) & ~isempty(num2str(Files.list{filenr}(tmp(1)+1)))
                    freqband = Files.list{filenr}(tmp(1)+1:end);
                    [~,~,~,freqband] = lab_filename(freqband);
                    tmp = strfind(freqband,'F');
                    if ~isempty(tmp) & length(tmp) > 1 & ~isempty(num2str(freqband(tmp(2)-1))) & ...
                            length(freqband) >= tmp(2)+1 & ~isempty(num2str(freqband(tmp(2)+1)))
                        tmp = strfind(freqband,'_');
                        if ~isempty(tmp)
                            freqband = freqband(1:tmp(1)-1);
                        end
                    else
                        tmp = strfind(freqband,'_');
                        if ~isempty(tmp) & ~isempty(num2str(freqband(tmp(1)-1))) & ...
                                length(freqband) >= tmp(1)+1 & ~isempty(num2str(freqband(tmp(1)+1)))
                            freqband(tmp(1)) = 'F';
                            if length(tmp) > 1
                                freqband = freqband(1:tmp(2)-1);                     
                            end
                        else
                            freqband = 'NoFreq';
                        end
                    end
                else
                    freqband = 'NoFreq';
                end
            end
            if ~exist('xlsout','var') | ~isfield(xlsout,freqband)
                for i = 1:size(MAT.Result.results,1)
                    for j = 1:numchans
                        resultstmp{j,i} = [MAT.Result.results{i,1} '_ch' num2str(j)]; %#ok<AGROW>
                    end
                end
                resultstmp = resultstmp(:);
                xlsout.(freqband) = cat(1,cellstr(''),resultstmp);
                clearvars resultstmp i j
                xlsoutmean.(freqband) = cat(1,cellstr('C1 R0'),MAT.Result.resultsmean);
                matrix.(freqband) = [];
            end
            if isempty(cfg.CollectGraph.subjectname) & isfield(MAT.Result,'patient')
                patient = MAT.Result.patient;
            elseif ~isempty(cfg.CollectGraph.subjectname)
                [patient,cfg.CollectConnect,skipprocessing] = lab_subjectname(Files.list{1,filenr},cfg.CollectGraph);
                if skipprocessing == 1
                    return
                end
            else
                patient = lab_subjectname(Files.list{1,filenr});
            end
            tmp = regexp(patient,'\d');
            if length(tmp) == length(patient)
                patient = ['P_' patient]; %#ok<AGROW>
            end
            if isempty(MAT.Result.data) & isfield(MAT.Result,'dataAVG') & ~isempty(MAT.Result.dataAVG)
                averagemode = 3;
            elseif ~isempty(MAT.Result.data) & (~isfield(MAT.Result,'dataAVG') | isempty(MAT.Result.dataAVG))
                if cfg.CollectGraph.mode < 3
                    averagemode = cfg.CollectGraph.mode;
                else
                    averagemode = 2;
                end
            elseif ~isempty(MAT.Result.dataAVG) & ~isempty(MAT.Result.data)
                averagemode = cfg.CollectGraph.mode;
            else
                skipfile = 1;
            end
        end
        if skipfile == 0
            if averagemode == 1
                disp('   collect results from single matrices')
                if size(MAT.Result.data,2) > 1
                    headertmp = cell(1,size(MAT.Result.data,2));
                    for i = 1:size(MAT.Result.data,2)
                        headertmp{1,i} = [patient '_' num2str(i)];
                    end
                else
                    headertmp{1,1} = patient;
                end
                datatmp = permute(MAT.Result.data,[1 3 2]);
                datatmp = reshape(datatmp,size(datatmp,1)*size(datatmp,2),size(datatmp,3));
                datatmpmean = MAT.Result.datamean;
                xlsout.(freqband){1,1} = ['C' num2str(size(MAT.Result.data,1)) ' R0'];
            elseif averagemode == 2
                headertmp = patient;
                datatmp = permute(MAT.Result.data,[1 3 2]);
                datatmp = reshape(datatmp,size(datatmp,1)*size(datatmp,2),size(datatmp,3));
                if isfield(cfg.CollectGraph,'doaveragenum') & ~isempty(cfg.CollectGraph.doaveragenum) & ...
                        size(datatmp,2) > cfg.CollectGraph.doaveragenum
                    disp(['   take average result of ' num2str(cfg.CollectGraph.doaveragenum) ' matrices'])
                    datatmp = mean(datatmp(:,1:cfg.CollectGraph.doaveragenum),2);
                    datatmpmean = mean(MAT.Result.datamean(:,1:cfg.CollectGraph.doaveragenum),2);
                else
                    disp(['   take average result of ' num2str(size(datatmp,2)) ' matrices'])
                    datatmp = mean(datatmp,2);
                    datatmpmean = mean(MAT.Result.datamean,2);
                end
                xlsout.(freqband){1,1} = ['C' num2str(size(MAT.Result.data,1)) ' R0'];
            elseif averagemode == 3
                headertmp = patient;
                disp('   collect results from average matrix')
                datatmp = permute(MAT.Result.dataAVG,[1 3 2]);
                datatmp = reshape(datatmp,size(datatmp,1)*size(datatmp,2),1);
                datatmpmean = MAT.Result.datameanAVG;
                xlsout.(freqband){1,1} = ['C' num2str(size(MAT.Result.dataAVG,1)) ' R0'];
            else
                skipfile = 1;
            end
        end
        if skipfile == 0
            if isfield(MAT.Result,'matrix') & ~isempty(MAT.Result.matrix)
                matrix.(freqband) = cat(3,matrix.(freqband),MAT.Result.matrix);
            end
            data = cat(1,headertmp,num2cell(datatmp));
            datamean = cat(1,headertmp,num2cell(datatmpmean));
            xlsout.(freqband) = [xlsout.(freqband) data];
            xlsoutmean.(freqband) = [xlsoutmean.(freqband) datamean];
            clearvars headertmp datatmp data MAT.Result doaverage averagemode
        end
    end
        
    % write Files
    freqbands = fieldnames(xlsout);
    Sort = [];
    for i = 1:length(freqbands)
        if strcmp(freqbands{i}(1),'F')
            tmp = freqbands{i}(2:end);
            tmp2 = strfind(tmp,'F');
            if ~isempty(tmp2)
                Sort(end+1) = str2num(tmp(1:tmp2(1)-1)); %#ok<ST2NM,AGROW>
            end
        end
    end
    if length(Sort) == length(freqbands);
        [~,Sort] = sort(Sort);
        freqbands = freqbands(Sort);
    end
    xlsout2 = {};
    xlsoutmean2 = {};
    Vars = {};
    Subjects = {};
    for i = 1:length(freqbands)
        if isempty(Vars)
            Vars = xlsout.(freqbands{i})(2:end,1);
            Subjects = xlsout.(freqbands{i})(1,2:end);
            xlsout2 = xlsout.(freqbands{i});
            VarsMean = xlsoutmean.(freqbands{i})(2:end,1);
            SubjectsMean = xlsoutmean.(freqbands{i})(1,2:end);
            xlsoutmean2 = xlsoutmean.(freqbands{i});
            if ~strcmp(freqbands{i},'NoFreq') & length(freqbands) > 1
                for j = 2:size(xlsout2,1)
                    xlsout2{j,1} = [freqbands{i} '_' xlsout2{j,1}];
                end
                for j = 2:size(xlsoutmean2,1)
                    xlsoutmean2{j,1} = [freqbands{i} '_' xlsoutmean2{j,1}];
                end
            end
        else
            [~,I1,I2] = intersect(Vars,xlsout.(freqbands{i})(2:end,1),'stable');
            [~,S1,S2] = intersect(Subjects,xlsout.(freqbands{i})(1,2:end),'stable');
            I1 = [1;I1+1];
            I2 = [1;I2+1];
            S1 = [1;S1+1];
            S2 = [1;S2+1];
            if ~isempty(I1)
                xlstmp = xlsout.(freqbands{i})(I2,S2);
                if ~strcmp(freqbands{i},'NoFreq')
                    for j = 2:size(xlstmp,1)
                        xlstmp{j,1} = [freqbands{i} '_' xlstmp{j,1}];
                    end
                end
                xlsout2 = cat(3,xlsout2(I1,S1,:),xlstmp);
                Vars = Vars(I1(2:end)-1);
                Subjects = Subjects(S1(2:end)-1);
            else
                disp(['Skip ' freqbands{i} ' no matching subjects or variables']);
            end
            [~,I1,I2] = intersect(VarsMean,xlsoutmean.(freqbands{i})(2:end,1),'stable');
            [~,S1,S2] = intersect(SubjectsMean,xlsoutmean.(freqbands{i})(1,2:end),'stable');
            I1 = [1;I1+1];
            I2 = [1;I2+1];
            S1 = [1;S1+1];
            S2 = [1;S2+1];
            if ~isempty(I1)
                xlstmp = xlsoutmean.(freqbands{i})(I2,S2);
                if ~strcmp(freqbands{i},'NoFreq')
                    for j = 2:size(xlstmp,1)
                        xlstmp{j,1} = [freqbands{i} '_' xlstmp{j,1}];
                    end
                end
                xlsoutmean2 = cat(3,xlsoutmean2(I1,S1,:),xlstmp);
                VarsMean = VarsMean(I1(2:end)-1);
                SubjectsMean = SubjectsMean(S1(2:end)-1);
            else
                disp(['Skip ' freqbands{i} ' no matching subjects or variables']);
            end
        end
    end
    
    Subjectstmp = xlsout2(1,:,1);
    xlsout2 = xlsout2(2:end,:,:);
    NumVars = size(xlsout2,1);
    xlsout2 = permute(xlsout2,[1 3 2]);
    xlsout2 = reshape(xlsout2,size(xlsout2,1)*size(xlsout2,2),size(xlsout2,3));
    xlsout2 = cat(1,Subjectstmp,xlsout2);
    if ~isempty(xlsout2{1,1})
        xlsout2{1,1} = [xlsout2{1,1} ' V' num2str(floor(NumVars/numchans))];
    else
        xlsout2{1,1} = ['C' num2str(numchans) ' R0 V' num2str(floor(NumVars/numchans))];
    end
    
    Subjectstmp = xlsoutmean2(1,:,1);
    xlsoutmean2 = xlsoutmean2(2:end,:,:);
    NumVars2 = size(xlsoutmean2,1);
    xlsoutmean2 = permute(xlsoutmean2,[1 3 2]);
    xlsoutmean2 = reshape(xlsoutmean2,size(xlsoutmean2,1)*size(xlsoutmean2,2),size(xlsoutmean2,3));
    xlsoutmean2 = cat(1,Subjectstmp,xlsoutmean2);
    xlsoutmean2{1,1} = ['C' num2str(NumVars2) ' R0'];
    
    %correct MST_ to MST-
    for i= 2:size(xlsout2,1)
        xlsout2{i,1} = regexprep(xlsout2{i,1},'MST_','MST-');
    end
    for i= 2:size(xlsoutmean2,1)
        xlsoutmean2{i,1} = regexprep(xlsoutmean2{i,1},'MST_','MST-');
    end
    
    % correct _Norm to -Norm
    for i= 2:size(xlsout2,1)
        xlsout2{i,1} = regexprep(xlsout2{i,1},'_Norm','-Norm');
    end
    for i= 2:size(xlsoutmean2,1)
        xlsoutmean2{i,1} = regexprep(xlsoutmean2{i,1},'_Norm','-Norm');
    end
    
    lab_write_xls(fullfile(Output_path,[Files.name '.xlsx']),xlsout2);
    lab_write_xls(fullfile(Output_path,[Files.name '_mean.xlsx']),xlsoutmean2);
        
    % Calculate MST
    if exist('matrix','var') & ~isempty(matrix)
        freqbands = fieldnames(matrix);
        for i = 1:length(freqbands)
            disp(['   Calculate MST for ' Files.name])
            matrix2 = mean(matrix.(freqbands{i}),3);
            MST = lab_MST(matrix2);
            if strcmp(freqbands{i},'NoFreq')
                Name = Files.name;
            else
                Name = [Files.name '_' freqbands{i}];
            end
            matrixfileout = fullfile(Output_path,[Name '_matrix.txt']);
            if exist(matrixfileout,'file')
                delete(matrixfileout);
            end
            dlmwrite(matrixfileout,matrix2,'delimiter','\t','precision', 6);
            
            matrixfileout = fullfile(Output_path,[Name '_MST.txt']);
            if exist(matrixfileout,'file')
                delete(matrixfileout);
            end
            dlmwrite(matrixfileout,MST,'delimiter','\t','precision', 6);
        end
    end
    
    % clear variables for next run
    clearvars matrix MST xlsout matrixfleout
end

% turn log-file off
diary off

return