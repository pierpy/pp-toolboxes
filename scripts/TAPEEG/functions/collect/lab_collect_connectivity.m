% Collect results from graphanalysis analysis
% Output is xls-files
%
% lab_collect_connectivitydata(cfg)
%
% written by F. Hatz 2013

function lab_collect_connectivity(cfg)
disp('Collect connectivity data of all files')

skipprocessing = 0;
if ~exist('cfg','var')
    cfg = [];
    skipselection = false;
else
    skipselection = true;
end
FilesAll = [];

if ~isfield(cfg,'CollectConnect') | ~isfield(cfg.CollectConnect,'searchfolder')
    [cfg,FilesAll,skipprocessing] = lab_set_collect_connectivity(cfg);
    if skipprocessing == 1
        return
    end
end

% turn log-file on
if exist(cfg.CollectConnect.searchfolder,'dir')
    diary(fullfile(cfg.CollectConnect.searchfolder,'CollectConnectivity.log'));
end

% search files
if isempty(FilesAll) | ~isfield(FilesAll,'list')
    FilesAll = lab_collect_connectivity_search(cfg,skipselection);
end
if isempty(FilesAll) | skipprocessing == 1
    disp('No connectivity results found')
    diary off
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
    [Filefirst,Filelast] = lab_find_lastperfolder(Files.list);
    R = [];
    for filenr = 1:size(Files.list,2)
        disp(['Collect connectivity data: ' lab_filename(Files.list{1,filenr})])
        load(Files.list{1,filenr});
        if filenr == size(Files.list,2)
            cfg.lastfile = true;
        else
            cfg.lastfile = false;
        end
        if filenr == 1
            cfg.firstfile = true;
        else
            cfg.firstfile = false;
        end
        cfg.lastfilefolder = Filelast(filenr);
        cfg.firstfilefolder = Filefirst(filenr);
        
        if exist('result','var')
            if isfield(result,'channels')
                header.channels = result.channels;
            else
                header = [];
            end
            if ~exist('patient','var') | ~isempty(cfg.CollectConnect.subjectname)
                [patient,cfg.CollectConnect,skipprocessing] = lab_subjectname(Files.list{1,filenr},cfg.CollectConnect);
                if skipprocessing == 1
                    return
                end
            end
            tmp = regexp(patient,'\d');
            if length(tmp) == length(patient)
                patient = ['P_' patient]; %#ok<AGROW>
            end
            clearvars tmp
            if ~isfield(R,'patient')
                R.patient = cellstr(patient);
            else
                R.patient = [R.patient cellstr(patient)];
            end
            cfg.Output_file = patient;
            if length(FilesAll) > 1 & isfield(Files,'name')
                cfg.Output_fileT = Files.name(6:end);
            else
                cfg.Output_fileT = '';
            end
            cfg.Output_fileS = cfg.Output_file;
            
            if isfield(cfg.CollectConnect,'randphase') & cfg.CollectConnect.randphase == true
                if isfield(Files,'listrand') & ~isempty(Files.listrand{1,filenr}) & exist(Files.listrand{1,filenr},'file')
                    tmp = load(Files.listrand{1,filenr});
                    if isfield(tmp,'result')
                        if ~isfield(cfg.CollectConnect,'methodrandphase') | isempty(cfg.CollectConnect.methodrandphase)
                            cfg.CollectConnect.methodrandphase = 'Threshold';
                        end
                        resultR = tmp.result;
                        result = correct_randphase(result,resultR,cfg.CollectConnect.valuerandphase,cfg.CollectConnect.value2randphase,cfg.CollectConnect.methodrandphase);
                    else
                        disp('   Skip correction by random phase, missing data')
                    end
                else
                    disp('   Skip correction by random phase, missing data')
                end
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
            if isfield(cfg.CollectConnect,'SelectVars') & cfg.CollectConnect.SelectVars == true
                if ~exist('SelectVars','var')
                    selection = listdlg('PromptString','Measures:','SelectionMode','multiple', ...
                        'ListString',Vars,'InitialValue',1:size(Vars,1),'CancelString','None','ListSize',[70 130]);
                    if isempty(selection)
                        return
                    else
                        pause(0.2);
                        SelectVars = Vars(selection,1);
                    end
                    clearvars selection
                end
                Vars = intersect(SelectVars,Vars);
            end
            for V = 1:length(Vars)
                if cfg.CollectConnect.doaverage == true
                    
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
                    Nmx = size(matrix,3);
                    if isfield(cfg.CollectConnect,'nummatrices') & ~isempty(cfg.CollectConnect.nummatrices) & ...
                            size(matrix,3) > cfg.CollectConnect.nummatrices
                        disp(['   Average matrices ' Vars{V} ' (' num2str(cfg.CollectConnect.nummatrices) ' matrices)'])
                        matrix = mean(matrix(:,:,1:cfg.CollectConnect.nummatrices),3);
                        Nmx = cfg.CollectConnect.nummatrices;
                    else
                        disp(['   Average matrices ' Vars{V} ' (' num2str(size(matrix,3)) ' matrices)'])
                        matrix = mean(matrix,3);
                    end
                    if isfield(cfg.CollectConnect,'MATRIX') & ~isempty(cfg.CollectConnect.MATRIX)
                        cfg.CollectConnect.Output_file = cfg.Output_file;
                        cfg.CollectConnect.Output_fileS = cfg.Output_fileS;
                        cfg.CollectConnect.Output_band = cfg.Output_fileT;
                        cfg.CollectConnect.patient = patient;
                        if length(Vars) > 1
                            warning off %#ok<WNOFF>
                            mkdir(fullfile(cfg.Output_filepath,Vars{V}));
                            warning on %#ok<WNON>
                            cfg.CollectConnect.Output_filepath = fullfile(cfg.Output_filepath,Vars{V});
                        else
                            cfg.CollectConnect.Output_filepath = cfg.Output_filepath;
                        end
                        if filenr == size(Files.list,2)
                            cfg.CollectConnect.lastfile = true;
                        else
                            cfg.CollectConnect.lastfile = false;
                        end
                        if filenr == 1
                            cfg.CollectConnect.firstfile = true;
                        else
                            cfg.CollectConnect.firstfile = false;
                        end
                        cfg.CollectConnect.lastfilefolder = Filelast(filenr);
                        cfg.CollectConnect.firstfilefolder = Filefirst(filenr);
                        header.measure = Vars{V};
                        [matrix,Rtmp,cfg.CollectConnect] = lab_process_matrix(matrix,header,cfg.CollectConnect,true,true);
                        if cfg.CollectConnect.MATRIX.grandaverage == 1
                            patient = cellstr(['GrandAverage_' Vars{V}]);
                            if isempty(cfg.Output_fileT)
                                cfg.Output_file = ['GrandAverage_' Vars{V}];
                            else  
                                cfg.Output_file = ['GrandAverage_' Vars{V} '_' cfg.Output_fileT];
                            end
                        end
                        if isfield(Rtmp,'channels')
                            header.channels = Rtmp.channels;
                        end
                    end
                    if ~isempty(matrix)
                        if isfield(cfg.CollectConnect,'GRAPH') & ~isempty(cfg.CollectConnect.GRAPH)
                            Output_file = cfg.Output_file;
                            Output_fileS = cfg.Output_fileS;
                            if isempty(cfg.Output_fileT)
                                cfg.Output_file = ['GraphResults_' Vars{V}];
                                cfg.Output_fileS = ['GraphResults_' Vars{V}];
                            else
                                cfg.Output_file = ['GraphResults_' Vars{V} '_' cfg.Output_fileT];
                                cfg.Output_fileS = ['GraphResults_' Vars{V} '_' cfg.Output_fileT];
                            end
                            Graph_folder = cfg.CollectConnect.GRAPH.folder;
                            if ~isempty(cfg.Output_fileT)
                                cfg.CollectConnect.GRAPH.folder = [cfg.CollectConnect.GRAPH.folder '_' cfg.Output_fileT '_' Vars{V}];
                            else
                                cfg.CollectConnect.GRAPH.folder = [cfg.CollectConnect.GRAPH.folder '_' Vars{V}];
                            end
                            [Rtmp,cfg.CollectConnect] = lab_graphanalysis(matrix,header,cfg.CollectConnect,cfg,cellstr(Output_file));
                            R.graph.(Vars{V}) = Rtmp;
                            cfg.CollectConnect.GRAPH.folder = Graph_folder;
                            cfg.Output_file = Output_file;
                            cfg.Output_fileS = Output_fileS;
                            clearvars Output_file Output_fileS Graph_folder
                        end
                        if ~isfield(R,'subjects') | ~isfield(R.subjects,Vars{V})
                            R.subjects.(Vars{V}) = cellstr(patient);
                        else
                            R.subjects.(Vars{V}) = [R.subjects.(Vars{V}) cellstr(patient)];
                        end
                        if ~isfield(R,'matrix') | ~isfield(R.matrix,Vars{V})
                            R.matrix.(Vars{V}) = matrix;
                            R.number.(Vars{V}) = Nmx;
                        else
                            R.matrix.(Vars{V}) = cat(3,R.matrix.(Vars{V}),matrix);
                            R.number.(Vars{V}) = cat(2,R.number.(Vars{V}),Nmx);
                        end
                        matrix(1:size(matrix,2)+1:end) = 0;
                        if ~isempty(setdiff(unique(matrix(:)),[0,1]))
                            degrees = sum(matrix,2) ./ (size(matrix,2)-1);
                        else
                            degrees = sum(matrix,2);
                        end
                        if ~isfield(R,'degrees') | ~isfield(R.degrees,Vars{V})
                            R.degrees.(Vars{V}) = degrees;
                        else
                            R.degrees.(Vars{V}) = cat(2,R.degrees.(Vars{V}),degrees);
                        end
                        if ~isfield(R,'channels') | ~isfield(R.channels,Vars{V})
                            if isfield(header,'channels') & size(matrix,1) == size(header.channels,1)
                                R.channels.(Vars{V}) = cellstr(header.channels);
                            else
                                R.channels.(Vars{V}) = cell(size(matrix,1),1);
                                for j = 1:size(matrix,2)
                                    R.channels.(Vars{V}){j,1} = ['Ch' num2str(j)];
                                end
                            end
                        end
                        clearvars matrix degrees
                    end
                else
                    matrix = result.(Vars{V});
                    Nmx = size(matrix,3);
                    if isfield(cfg.CollectConnect,'nummatrices') & ~isempty(cfg.CollectConnect.nummatrices) & ...
                            size(matrix,3) > cfg.CollectConnect.nummatrices
                        matrix = matrix(:,:,1:cfg.CollectConnect.nummatrices);
                        Nmx = cfg.CollectConnect.nummatrices;
                    end
                    if isfield(cfg.CollectConnect,'binary') & cfg.CollectConnect.binary == true
                        if isempty(cfg.CollectConnect.binarythreshold)
                            cfg.CollectConnect.binarythreshold = 0.5;
                        end
                        if isempty(cfg.CollectConnect.binarymode)
                            cfg.CollectConnect.binarymode = 'fixed';
                        end
                        matrix = lab_matrix2binary(matrix,cfg.CollectConnect.binarythreshold,cfg.CollectConnect.binarymode);
                    end
                    patienttmp = {};
                    for i = 1:Nmx
                        patienttmp{1,i} = [patient '_' num2str(i)]; %#ok<AGROW>
                    end
                    if isfield(cfg.CollectConnect,'MATRIX') & ~isempty(cfg.CollectConnect.MATRIX)
                        % cfg.CollectConnect.Output_file = cfg.Output_file;
                        % cfg.CollectConnect.Output_fileS = cfg.Output_file;
                        % cfg.CollectConnect.Output_band = cfg.Output_fileT;
                        % if length(Vars) > 1
                        %     warning off %#ok<WNOFF>
                        %     mkdir(fullfile(cfg.Output_filepath,Vars{V}));
                        %     warning on %#ok<WNON>
                        %     cfg.CollectConnect.Output_filepath = fullfile(cfg.Output_filepath,Vars{V});
                        % else
                        %     cfg.CollectConnect.Output_filepath = cfg.Output_filepath;
                        % end
                        if filenr == size(Files.list,2)
                            cfg.CollectConnect.lastfile = true;
                        else
                            cfg.CollectConnect.lastfile = false;
                        end
                        if filenr == 1
                            cfg.CollectConnect.firstfile = true;
                        else
                            cfg.CollectConnect.firstfile = false;
                        end
                        cfg.CollectConnect.lastfilefolder = Filelast(filenr);
                        cfg.CollectConnect.firstfilefolder = Filefirst(filenr);
                        header.measure = Vars{V};
                        [matrix,Rtmp,cfg.CollectConnect] = lab_process_matrix(matrix,header,cfg.CollectConnect,true,true);
                        if cfg.CollectConnect.MATRIX.grandaverage == 1
                            patienttmp = cellstr(['GrandAverage_' Vars{V}]);
                            if isempty(cfg.Output_fileT)
                                cfg.Output_file = ['GrandAverage_' Vars{V}];
                            else  
                                cfg.Output_file = ['GrandAverage_' Vars{V} '_' cfg.Output_fileT];
                            end
                        end
                        if isfield(Rtmp,'channels')
                            header.channels = Rtmp.channels;
                        end
                    end
                    if ~isempty(matrix)
                        if isfield(cfg.CollectConnect,'GRAPH') & ~isempty(cfg.CollectConnect.GRAPH)
                            Output_file = cfg.Output_file;
                            Output_fileS = cfg.Output_fileS;
                            if isempty(cfg.Output_fileT)
                                cfg.Output_file = ['GraphResults_' Vars{V}];
                                cfg.Output_fileS = ['GraphResults_' Vars{V}];
                            else
                                cfg.Output_file = ['GraphResults_' Vars{V} '_' cfg.Output_fileT];
                                cfg.Output_fileS = ['GraphResults_' Vars{V} '_' cfg.Output_fileT];
                            end
                            Graph_folder = cfg.CollectConnect.GRAPH.folder;
                            if ~isempty(cfg.Output_fileT)
                                cfg.CollectConnect.GRAPH.folder = [cfg.CollectConnect.GRAPH.folder '_' cfg.Output_fileT '_' Vars{V}];
                            else
                                cfg.CollectConnect.GRAPH.folder = [cfg.CollectConnect.GRAPH.folder '_' Vars{V}];
                            end
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
                            cfg.CollectConnect.GRAPH.folder = Graph_folder;
                            cfg.Output_file = Output_file;
                            cfg.Output_fileS = Output_fileS;
                            clearvars Output_file Output_fileS Output_file2 Graph_folder
                        end
                        if ~isfield(R,'subjects') | ~isfield(R.subjects,Vars{V})
                            R.subjects.(Vars{V}) = patienttmp;
                        else
                            R.subjects.(Vars{V}) = [R.subjects.(Vars{V}) patienttmp];
                        end
                        if ~isfield(R,'matrix') | ~isfield(R.matrix,Vars{V})
                            R.matrix.(Vars{V}) = matrix;
                        else
                            R.matrix.(Vars{V}) = cat(3,R.matrix.(Vars{V}),matrix);
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
                        if ~isfield(R,'degrees') | ~isfield(R.degrees,Vars{V})
                            R.degrees.(Vars{V}) = degrees;
                        else
                            R.degrees.(Vars{V}) = cat(2,R.degrees.(Vars{V}),degrees);
                        end
                        if ~isfield(R,'channels') | ~isfield(R.channels,Vars{V})
                            if isfield(header,'channels') & size(matrix,1) == size(header.channels,1)
                                R.channels.(Vars{V}) = cellstr(header.channels);
                            else
                                R.channels.(Vars{V}) = cell(size(matrix,1),1);
                                for j = 1:size(matrix,2)
                                    R.channels.(Vars{V}){j,1} = ['Ch' num2str(j)];
                                end
                            end
                        end
                        clearvars matrix degrees patienttmp
                    end
                end
            end
        end
    end
    
    % Store Result
    if isfield(Files,'freqband')
        Result.(Files.freqband) = R;
    elseif isfield(Files,'name')
        Result.(Files.name) = R;
    else
        Result.('NoFreq') = R;
    end
    
    if isfield(cfg.CollectConnect,'WriteMatrices') & cfg.CollectConnect.WriteMatrices == true
        disp('Write Matrices')
        Vars = fieldnames(R.matrix);
        for V = 1:length(Vars)
            if ~isempty(R.matrix.(Vars{V})) & isnumeric(R.matrix.(Vars{V})) & ...
                    size(R.matrix.(Vars{V}),1) == size(R.matrix.(Vars{V}),2)
                if isfield(Files,'name')
                    if length(Vars) > 1
                        Output_filepath2 = fullfile(cfg.Output_filepath,[Files.name '_' Vars{V}]);
                    else
                        Output_filepath2 = fullfile(cfg.Output_filepath,Files.name);
                    end
                else
                    if length(Vars) > 1
                        Output_filepath2 = fullfile(cfg.Output_filepath,Vars{V});
                    else
                        Output_filepath2 = cfg.Output_filepath;
                    end
                end
                warning off %#ok<WNOFF>
                mkdir(Output_filepath2);
                warning on %#ok<WNON>
                for i = 1:size(R.matrix.(Vars{V}),3)
                    matrixfileout = fullfile(Output_filepath2,[R.subjects.(Vars{V}){i} '_' Vars{V} '_matrix.txt']);
                    if exist(matrixfileout,'file')
                        delete(matrixfileout);
                    end
                    dlmwrite(matrixfileout,R.matrix.(Vars{V})(:,:,i),'delimiter','\t','precision', 6);
                end
            end
        end
    end
    
    if isfield(cfg.CollectConnect,'PLOT') & ~isempty(cfg.CollectConnect.PLOT)
        Vars = fieldnames(R.matrix);
        for V = 1:length(Vars)
            if ~isempty(R.matrix.(Vars{V})) & isnumeric(R.matrix.(Vars{V})) & ...
                    size(R.matrix.(Vars{V}),1) == size(R.matrix.(Vars{V}),2)
                if isfield(Files,'name')
                    if length(Vars) > 1
                        Output_filepath2 = fullfile(cfg.Output_filepath,[Files.name '_' Vars{V}]);
                    else
                        Output_filepath2 = fullfile(cfg.Output_filepath,Files.name);
                    end
                else
                    if length(Vars) > 1
                        Output_filepath2 = fullfile(cfg.Output_filepath,Vars{V});
                    else
                        Output_filepath2 = cfg.Output_filepath;
                    end
                end
                warning off %#ok<WNOFF>
                mkdir(Output_filepath2);
                warning on %#ok<WNON>
                for i = 1:size(R.matrix.(Vars{V}),3)
                    matrixfileout = fullfile(Output_filepath2,[R.subjects.(Vars{V}){i} '_' Vars{V} '_matrix.tif']);
                    cfg.CollectConnect.PLOT = plot_matrix(R.matrix.(Vars{V})(:,:,i),matrixfileout,cfg.CollectConnect.PLOT);
                end
            end
        end
    end
end
if exist('Result','var')
    save(fullfile(cfg.Output_filepath,'ResultsConnectivity.mat'),'Result');
end

freqbands = fieldnames(Result);
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

% Write Nr-Files
if isfield(Result.(freqbands{1}),'number')
    disp('Write Nr-Files')
    Vars = fieldnames(Result.(freqbands{1}).number);
    for V = 1:length(Vars)
        fileout = fullfile(cfg.Output_filepath,[upper(Vars{V}) '_Nr.xlsx']);
        if exist(fileout,'file')
            delete(fileout);
        end
        xlsout = {};
        subjects = {};
        measures = {};
        for F = 1:length(freqbands)
            if isfield(Result.(freqbands{F}),'number') & isfield(Result.(freqbands{F}).number,Vars{V}) & ...
                    ~isempty(Result.(freqbands{F}).number.(Vars{V}))
                if isempty(subjects)
                    subjects = Result.(freqbands{1}).subjects.(Vars{V});
                    measures = freqbands(F);
                    xlsout = num2cell(Result.(freqbands{1}).number.(Vars{V}));
                else
                    [subjects,S1,S2] = intersect(subjects,Result.(freqbands{1}).subjects.(Vars{V}),'stable');
                    xlsout = cat(1,xlsout(:,S1),num2cell(Result.(freqbands{1}).number.(Vars{V})(:,S2)));
                    measures = cat(1,measures,freqbands(F));
                end
            end
        end
        xlsout = cat(1,subjects,xlsout);
        xlsout = cat(2,cat(1,{'C1 R0'},measures),xlsout);
        lab_write_xls(fileout,xlsout');
        clearvars xlsout
    end
    clearvars Vars V F
end

% Write XLS-Files with degrees
if isfield(cfg.CollectConnect,'WriteDegrees') & cfg.CollectConnect.WriteDegrees == true & isfield(Result.(freqbands{1}),'degrees')
    disp('Write Degrees')
    Vars = fieldnames(Result.(freqbands{1}).degrees);
    for V = 1:length(Vars)
        fileout = fullfile(cfg.Output_filepath,[upper(Vars{V}) '_Degree.xlsx']);
        if exist(fileout,'file')
            delete(fileout);
        end
        xlsout = {};
        subjects = {};
        measures = {};
        ICluster = [];
        for F = 1:length(freqbands)
            if isfield(Result.(freqbands{F}),'degrees') & isfield(Result.(freqbands{F}).degrees,Vars{V}) & ...
                    ~isempty(Result.(freqbands{F}).degrees.(Vars{V}))
                Mtmp = Result.(freqbands{F}).channels.(Vars{V})(:);
                for i = 1:length(Result.(freqbands{F}).channels.(Vars{V}))
                    Mtmp{i} = [freqbands{F} '_' Mtmp{i}];
                end
                if isempty(subjects)
                    subjects = Result.(freqbands{1}).subjects.(Vars{V});
                    measures = Mtmp;
                    xlsout = num2cell(Result.(freqbands{1}).degrees.(Vars{V}));
                    ICluster = length(Mtmp);
                else
                    [subjects,S1,S2] = intersect(subjects,Result.(freqbands{1}).subjects.(Vars{V}),'stable');
                    xlsout = cat(1,xlsout(:,S1),num2cell(Result.(freqbands{1}).degrees.(Vars{V})(:,S2)));
                    measures = cat(1,measures,Mtmp);
                    ICluster = [ICluster length(Mtmp)]; %#ok<AGROW>
                end
            end
        end
        xlsout = cat(1,subjects,xlsout);
        if min(ICluster) == max(ICluster)
            xlsout = cat(2,cat(1,{['C' num2str(min(ICluster)) ' R0']},measures),xlsout);
        else
            xlsout = cat(2,cat(1,{'C1 R0'},measures),xlsout);
        end
        lab_write_xls(fileout,xlsout);
        clearvars xlsout
    end
    clearvars Vars V F
end
    
% Write XLS-Files with connections
if isfield(cfg.CollectConnect,'WriteConnections') & cfg.CollectConnect.WriteConnections == true
    disp('Write Connections')
    Vars = fieldnames(Result.(freqbands{1}).matrix);
    for V = 1:length(Vars)
        fileout = fullfile(cfg.Output_filepath,[upper(Vars{V}) '_Connections.xlsx']);
        if exist(fileout,'file')
            delete(fileout);
        end
        xlsout = {};
        subjects = {};
        measures = {};
        ICluster = [];
        for F = 1:length(freqbands)
            if isfield(Result.(freqbands{F}),'matrix') & isfield(Result.(freqbands{F}).matrix,Vars{V}) & ...
                    ~isempty(Result.(freqbands{F}).matrix.(Vars{V}))
                [connections,channels] = lab_extract_tril(Result.(freqbands{F}).matrix.(Vars{V}),Result.(freqbands{F}).channels.(Vars{V}));
                for i = 1:length(channels)
                    channels{i} = [freqbands{F} '_' channels{i}];
                end
                if isempty(subjects)
                    subjects = Result.(freqbands{1}).subjects.(Vars{V});
                    measures = channels;
                    xlsout = num2cell(connections);
                    ICluster = length(channels);
                else
                    [subjects,S1,S2] = intersect(subjects,Result.(freqbands{1}).subjects.(Vars{V}),'stable');
                    xlsout = cat(1,xlsout(:,S1),num2cell(connections(:,S2)));
                    measures = cat(1,measures,channels);
                    ICluster = [ICluster length(channels)]; %#ok<AGROW>
                end
            end
        end
        xlsout = cat(1,subjects,xlsout);
        if min(ICluster) == max(ICluster)
            xlsout = cat(2,cat(1,{['C' num2str(min(ICluster)) ' R0']},measures),xlsout);
        else
            xlsout = cat(2,cat(1,{'C1 R0'},measures),xlsout);
        end
        lab_write_xls(fileout,xlsout);
        clearvars xlsout
    end
    clearvars Vars V F
end

if isfield(cfg.CollectConnect,'Kmeans') & ~isempty(cfg.CollectConnect.Kmeans)
    Output_filepath2 = fullfile(cfg.Output_filepath,'Kmeans');
    warning off %#ok<WNOFF>
    mkdir(Output_filepath2);
    warning on %#ok<WNON>
    cfg.CollectConnect.Kmeans.searchfolder = Output_filepath2;
    lab_matrix_clustering(Result,cfg.CollectConnect.Kmeans);
end

% turn log-file off
diary off

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
            Mtmp = mean(repmat(Mtmp,[1 1 Mrun]) >= Rtmp,3);
        end
        Mtmp(1:size(Mtmp,1)+1:end) = 0;
        Moutput(:,:,m) = Mtmp;
    end
end

function settings = plot_matrix(matrix,matrix_file,settings)
    disp(['    Plot matrix: ' matrix_file])
    settings.handleF = figure('Visible','off','Color',[1 1 1],'Renderer','zbuffer','Menubar','none');
    settings.PLOT_file = [];
    settings.close = 0;
    PLOT.ColorE = lab_create_cmap(settings.ColorE);
    PLOT.SizeE = settings.SizeE;
    PLOT.Color = lab_create_cmap(settings.Color);
    PLOT.Size = settings.Size;
    if isfield(settings,'nodemethod')
        nodes = lab_plot_create_nodes(matrix,settings);
    else
        nodes = zeros(1,size(matrix,1));
    end
    matrix(1:size(matrix,1)+1:end) = nodes;
    if max(nodes(:)) == 0 & min(nodes(:)) == 0
        PLOT.MinValue = 0;
        PLOT.MaxValue = 1;
    else
        PLOT.MinValue = min(nodes(:));
        PLOT.MaxValue = max(nodes(:));
    end
    settings = lab_plot_chans(matrix,PLOT,settings);
    lab_print_figure(matrix_file,settings.handleF);
    close(settings.handleF);
end

