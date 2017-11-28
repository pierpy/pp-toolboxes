% Collect results from graphanalysis analysis
% Output is xls-files
%
% lab_collect_connectivitydata_random(cfg)
%
%    modified lab_collect_connectivity to collect data from random-phase
%    only
%
% written by F. Hatz 2013

function lab_collect_connectivity_random(cfg)
disp('Collect connectivity data of all files')

if ~exist('cfg','var')
    cfg = [];
end

if ~isfield(cfg,'CollectConnect') | ~isfield(cfg.CollectConnect,'searchfolder') | ~exist(cfg.CollectConnect.searchfolder,'dir')
    [cfg,skipprocessing] = lab_set_collect_connectivity(cfg);
    if skipprocessing == 1
        return
    end
end

% search files
disp('Search for result files')
[cfg,FilesAll,skipprocessing] = lab_collect_connectivity_search_random(cfg);
if isempty(FilesAll) | skipprocessing == 1
    disp('no connectivity results found')
    return
end

% Create Output-folder
warning off %#ok<WNOFF>
if ~isempty(cfg.CollectConnect.outputfolder)
    mkdir(fullfile(cfg.CollectConnect.searchfolder,cfg.CollectConnect.outputfolder));
    cfg.Output_filepath = fullfile(cfg.CollectConnect.searchfolder,cfg.CollectConnect.outputfolder);
else
    mkdir(fullfile(cfg.CollectConnect.searchfolder,'ConnectivityAnalysis'));
    cfg.Output_filepath = fullfile(cfg.CollectConnect.searchfolder,'ConnectivityAnalysis');
end
warning on %#ok<WNON>

for connectnr = 1:size(FilesAll,2)
    Files = FilesAll(1,connectnr);
    R = [];
    xlsout = [];
    for filenr = 1:size(Files.list,2)
        disp(['Collect connectivity data: ' lab_filename(Files.list{1,filenr})])
        load(Files.list{1,filenr});
        if filenr == size(Files.list,2)
            cfg.lastfile = true;
        else
            cfg.lastfile = false;
        end
        
        if exist('result','var')
            if isfield(result,'channels')
                header.channels = result.channels;
            else
                header = [];
            end
            if ~exist('patient','var') | isempty(patient)
                patient = Files.list{1,filenr};
            end
            [patient,cfg.CollectConnect,skipprocessing] = lab_subjectname(patient,cfg.CollectConnect);
            if skipprocessing == 1
                return
            end
            tmp = regexp(patient,'\d');
            if length(tmp) == length(patient)
                patient = ['P_' patient];
            end
            clearvars tmp
            if ~isfield(R,'patient')
                R.patient = cellstr(patient);
            else
                R.patient = [R.patient cellstr(patient)];
            end
            cfg.Output_file = patient;
            cfg.Output_fileS = patient;
            
            if isfield(cfg.CollectConnect,'randphase') & cfg.CollectConnect.randphase == true
                disp('   Skip correction by random phase, missing data')
            end
            
            Vars = fieldnames(result);
            i = 1;
            while i <= length(Vars)
                if isempty(result.(Vars{i})) | ~isnumeric(result.(Vars{i})) | size(result.(Vars{i}),1) ~= size(result.(Vars{i}),2)
                    if i == 1 & length(Vars) > 1
                        Vars = Vars(2:end,1);
                    elseif i < length(Vars)
                        Vars = cat(1,Vars(1:i-1,1),Vars(i+1:end,1));
                    else
                        Vars = Vars(1:i-1,1);
                    end
                else
                    i = i + 1;
                end
            end
            clearvars i
            for V = 1:length(Vars)
                result.(Vars{V}) = result.(Vars{V})(:,:,1:100:end);
                if cfg.CollectConnect.doaverage == true
                    disp('   Average matrices')
                    if isfield(cfg.CollectConnect,'binary') & cfg.CollectConnect.binary == true
                        if isempty(cfg.CollectConnect.binarythreshold)
                            cfg.CollectConnect.binarythreshold = 0.5;
                        end
                        if isempty(cfg.CollectConnect.binarymode)
                            cfg.CollectConnect.binarymode = 'fixed';
                        end
                        matrix = lab_matrix2binary(result.(Vars{V}),cfg.CollectConnect.binarythreshold,cfg.CollectConnect.binarymode);
                    else
                        matrix = result.(Vars{V});
                    end
                    matrix = mean(matrix,3);
                    if isfield(cfg.CollectConnect,'MATRIX') & ~isempty(cfg.CollectConnect.MATRIX)
                        cfg.CollectConnect.Output_file = cfg.Output_file;
                        cfg.CollectConnect.Output_fileS = cfg.Output_file;
                        cfg.CollectConnect.Output_filepath = cfg.Output_filepath;
                        if filenr == size(Files.list,2)
                            cfg.CollectConnect.lastfile = true;
                        else
                            cfg.CollectConnect.lastfile = false;
                        end
                        [matrix,Rtmp,cfg.CollectConnect] = lab_process_matrix(matrix,header,cfg.CollectConnect,true,true);
                        if cfg.CollectConnect.MATRIX.grandaverage == 1
                            patient = cellstr('GrandAverage');
                            cfg.Output_file = 'GrandAverage';
                        end
                        if isfield(Rtmp,'channels')
                            header.channels = Rtmp.channels;
                        end
                    end
                    if ~isempty(matrix)
                        if isfield(cfg.CollectConnect,'GRAPH') & ~isempty(cfg.CollectConnect.GRAPH)
                            Output_file = cfg.Output_file;
                            Output_fileS = cfg.Output_fileS;
                            cfg.Output_file = 'GraphResults';
                            cfg.Output_fileS = 'GraphResults';
                            [Rtmp,cfg.CollectConnect] = lab_graphanalysis(matrix,header,cfg.CollectConnect,cfg,cellstr(cfg.Output_file));
                            R.graph.(Vars{V}) = Rtmp;
                            cfg.Output_file = Output_file;
                            cfg.Output_fileS = Output_fileS;
                            clearvars Output_file Output_fileS
                        end
                        if ~isfield(R,'store') | ~isfield(R.store,Vars{V})
                            R.store.(Vars{V}) = cellstr(patient);
                        else
                            R.store.(Vars{V}) = [R.store.(Vars{V}) cellstr(patient)];
                        end
                        if ~isfield(R,Vars{V})
                            R.(Vars{V}) = matrix;
                        else
                            R.(Vars{V}) = cat(3,R.(Vars{V}),matrix);
                        end
                        matrix(1:size(matrix,2)+1:end) = 0;
                        if ~isempty(setdiff(unique(matrix(:)),[0,1]))
                            degrees = sum(matrix,2) ./ (size(matrix,2)-1);
                        else
                            degrees = sum(matrix,2);
                        end
                        degrees = cat(1,cellstr(patient),num2cell(degrees));
                        if ~isfield(xlsout,Vars{V})
                            for j = 1:size(matrix,2)
                                Var{j,1} = [upper(Vars{V}) '_ch' num2str(j)];
                            end
                            xlsout.(Vars{V}) = cat(1,{''},Var);
                        end
                        xlsout.(Vars{V}) = cat(2,xlsout.(Vars{V}),degrees);
                        clearvars degrees
                    end
                else
                    patienttmp = {};
                    for i = 1:size(result.(Vars{V}),3)
                        patienttmp{1,i} = [patient '_' num2str(i)];
                    end
                    if isfield(cfg.CollectConnect,'binary') & cfg.CollectConnect.binary == true
                        if isempty(cfg.CollectConnect.binarythreshold)
                            cfg.CollectConnect.binarythreshold = 0.5;
                        end
                        if isempty(cfg.CollectConnect.binarymode)
                            cfg.CollectConnect.binarymode = 'fixed';
                        end
                        matrix = lab_matrix2binary(result.(Vars{V}),cfg.CollectConnect.binarythreshold,cfg.CollectConnect.binarymode);
                    else
                        matrix = result.(Vars{V});
                    end
                    if isfield(cfg.CollectConnect,'MATRIX') & ~isempty(cfg.CollectConnect.MATRIX)
                        cfg.CollectConnect.Output_file = cfg.Output_file;
                        cfg.CollectConnect.Output_fileS = cfg.Output_file;
                        cfg.CollectConnect.Output_filepath = cfg.Output_filepath;
                        if filenr == size(Files.list,2)
                            cfg.CollectConnect.lastfile = true;
                        else
                            cfg.CollectConnect.lastfile = false;
                        end
                        [matrix,Rtmp,cfg.CollectConnect] = lab_process_matrix(matrix,header,cfg.CollectConnect,true,true);
                        if cfg.CollectConnect.MATRIX.grandaverage == 1
                            patienttmp = cellstr('GrandAverage');
                            cfg.Output_file = 'GrandAverage';
                        end
                        if isfield(Rtmp,'channels')
                            header.channels = Rtmp.channels;
                        end
                    end
                    if ~isempty(matrix)
                        if isfield(cfg.CollectConnect,'GRAPH') & ~isempty(cfg.CollectConnect.GRAPH)
                            Output_file = cfg.Output_file;
                            Output_fileS = cfg.Output_fileS;
                            cfg.Output_file = 'GraphResults';
                            cfg.Output_fileS = 'GraphResults';
                            if size(matrix,3) > 1
                                Output_file2 = {};
                                for Nmatrix = 1:size(matrix,3)
                                    Output_file2 = cat(2,Output_file2,cellstr([cfg.Output_file '_' num2str(Nmatrix)]));
                                end
                            else
                                Output_file2 = cellstr(cfg.Output_file);
                            end
                            [Rtmp,cfg.CollectConnect] = lab_graphanalysis(matrix,header,cfg.CollectConnect,cfg,Output_file2);
                            R.graph.(Vars{V}) = Rtmp;
                            cfg.Output_file = Output_file;
                            cfg.Output_fileS = Output_fileS;
                            clearvars Output_file Output_fileS Output_file2
                        end
                        if ~isfield(R,'store') | ~isfield(R.store,Vars{V})
                            R.store.(Vars{V}) = patienttmp;
                        else
                            R.store.(Vars{V}) = [R.store.(Vars{V}) patienttmp];
                        end
                        if ~isfield(R,Vars{V})
                            R.(Vars{V}) = matrix;
                        else
                            R.(Vars{V}) = cat(3,R.(Vars{V}),matrix);
                        end
                        for i = 1:size(matrix,3)
                            tmp = matrix(:,:,i);
                            tmp(1:size(tmp,1)+1:end) = 0;
                            matrix(:,:,i) = tmp;
                        end
                        clearvars i tmp
                        if ~isempty(setdiff(unique(matrix(:)),[0,1]))
                            degrees = sum(permute(matrix,[1 3 2]),3) ./ (size(matrix,1)-1);
                        else
                            degrees = sum(permute(matrix,[1 3 2]),3);
                        end
                        degrees = cat(1,patienttmp,num2cell(degrees));
                        if ~isfield(xlsout,Vars{V})
                            for j = 1:size(matrix,2)
                                Var{j,1} = [upper(Vars{V}) '_ch' num2str(j)];
                            end
                            xlsout.(Vars{V}) = cat(1,{''},Var);
                        end
                        xlsout.(Vars{V}) = cat(2,xlsout.(Vars{V}),degrees);
                        clearvars patienttmp degrees
                    end
                end
            end
        end
    end

    % Write Files
    Vars = fieldnames(xlsout);
    for V = 1:length(Vars)
        if ~isempty(xlsout.(Vars{V})) & size(xlsout.(Vars{V}),2) > 1
            fileout = fullfile(cfg.Output_filepath,[Files.name '_' upper(Vars{V}) '.xlsx']);
            lab_write_xls(fileout,xlsout.(Vars{V}));
        end
    end
    warning off %#ok<WNOFF>
    mkdir(fullfile(cfg.Output_filepath,Files.name));
    warning on %#ok<WNON>
    Output_filepath2 = fullfile(cfg.Output_filepath,Files.name);
    Vars = fieldnames(R);
    for V = 1:length(Vars)
        if ~isempty(R.(Vars{V})) & isnumeric(R.(Vars{V})) & size(R.(Vars{V}),1) == size(R.(Vars{V}),2)
            for i = 1:size(R.(Vars{V}),3)
                matrixfileout = fullfile(Output_filepath2,[R.store.(Vars{V}){i} '_' Vars{V} '_matrix.txt']);
                if exist(matrixfileout,'file')
                    delete(matrixfileout);
                end
                dlmwrite(matrixfileout,R.(Vars{V})(:,:,i),'delimiter','\t','precision', 6);
            end
        end
    end
    Result.(Files.name) = R;
end
if exist('Result','var')
    save(fullfile(cfg.Output_filepath,'ResultsConnectivity.mat'),'Result');
end

end

function Result = correct_randphase(Result,ResultR,Value,Value2,Method)
   if ~strcmp(Method,'P-Value')
       disp(['   Correct by phase randomization: method: ' Method ' / threshold: ' num2str(Value) ' ' Value2])
   else
       disp(['   Correct by phase randomization: method: ' Method])
   end
   variables = fieldnames(Result);
   for i = 1:length(variables)
       if ~isstruct(Result.(variables{i}))
           if isnumeric(Result.(variables{i})) & size(Result.(variables{i}),1) > 1 & ...
                       size(Result.(variables{i}),2) > 1 & size(Result.(variables{i}),1) == size(Result.(variables{i}),2)
               Result.(variables{i}) = do_calc(Result.(variables{i}),ResultR.(variables{i}),Value,Value2,Method);
           end
       else
           variables2 = fieldnames(Result.(variables{i}));
           for j = 1:length(variables2)
               if isnumeric(Result.(variables{i}).(variables2{j})) & size(Result.(variables{i}).(variables2{j}),1) > 1 & ...
                       size(Result.(variables{i}).(variables2{j}),2) > 1 & ...
                       size(Result.(variables{i}).(variables2{j}),1) == size(Result.(variables{i}).(variables2{j}),2)
                   Result.(variables{i}).(variables2{j}) = do_calc(Result.(variables{i}).(variables2{j}), ...
                       ResultR.(variables{i}).(variables2{j}),Value,Value2,Method);
               end
           end
       end
   end
   
   
end

function Moutput = do_calc(Minput,Rinput,Value,Value2,Method)
    Mnum = size(Minput,3);
    Mrun = floor(size(Rinput,3) / Mnum);
    Moutput = zeros(size(Minput));
    if strcmp(Value2,'percent')
        if isempty(Value)
            disp('    Abort: skip rand phase, wrong settings');
            Moutput = Minput;
            return
        end
        if ceil((Value/100)*Mrun) > (Value/100)*Mrun
            Value = ceil((Value/100)*Mrun);
        else
            Value = ceil((Value/100)*Mrun) + 1;
        end
        if Value < 1;Value = 1;end
        if Value > Mrun;Value = Mrun;end
    end
    for m = 1:Mnum
        if strcmp(Value2,'std')
            Rtmp = mean(Rinput(:,:,(m-1)*Mrun+1:m*Mrun),3) + ...
                std(Rinput(:,:,(m-1)*Mrun+1:m*Mrun),[],3)*Value;
        elseif strcmp(Value2,'percent')
            Rtmp = sort(Rinput(:,:,(m-1)*Mrun+1:m*Mrun),3);
            Rtmp = Rtmp(:,:,Value);
        elseif strcmp(Method,'P-Value')
            Rtmp = Rinput(:,:,(m-1)*Mrun+1:m*Mrun);
        else
            disp('    Abort: skip rand phase, wrong settings');
            Moutput = Minput;
            return
        end
        Mtmp = Minput(:,:,m);
        if strcmp(Method,'Threshold')
            Mtmp((Mtmp-Rtmp)<0) = 0;
        elseif strcmp(Method,'Diff')
            Mtmp = Mtmp - Rtmp;
            Mtmp(Mtmp<0) = 0;
        else
            Mtmp = 1 - mean(repmat(Mtmp,[1 1 Mrun]) >= Rtmp,3);
        end
        Mtmp(1:size(Mtmp,1)+1:end) = 0;
        Moutput(:,:,m) = Mtmp;
    end
end