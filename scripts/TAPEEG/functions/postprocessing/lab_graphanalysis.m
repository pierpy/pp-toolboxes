% Graph analysis
%
% [Result,settings] = lab_graphanalysis(matrixAll,header,settings,timestamp)
%
% matrixAll = array (nodes x nodes x timepoints)
% header    = see lab_create_header
% settings  = structure with config (optional)
%             settings.Output_file = 'Result.xls';
%             settings.Output_filepath
% timestamp = label for every matrix (if multiple timepoints)
%
% written by F. Hatz 2012

function [Result,settings,cfg] = lab_graphanalysis(matrixAll,header,settings,cfg,timestamp,doIncrease)

disp('   Calculate graph measures: ')

if ~exist('doIncrease','var')
    doIncrease = false;
end
if ~exist('cfg','var')
    cfg = [];
end
if ~exist('settings','var')
    settings = [];
end
if ~exist('header','var')
    header = [];
end
if ~exist('matrixAll','var')
    matrixAll = lab_read_data;
    if ~isnumeric(matrixAll) | size(matrixAll,1) ~= size(matrixAll,2)
        disp('Abort: wrong input data')
        return
    end
end

if ~isfield(settings,'GRAPH') | ~isfield(settings.GRAPH,'measure')
    [settings,skipprocessing] = lab_set_graphanalysis(settings,size(matrixAll,3),cfg);
    if skipprocessing == 1
        Result = [];
        return
    end
end

if isfield(settings,'isdpli') & settings.isdpli == true
    matrixAll = (matrixAll - 0.5) * 2;
end

if exist('timestamp','var') & ~isempty(timestamp)
    filenames = timestamp(:)';
elseif size(matrixAll,3) == 1
    if isfield(cfg,'EEG_fileS')
        filenames = cellstr(cfg.EEG_fileS);
    elseif isfield(cfg,'Output_fileS')
        filenames = cellstr(cfg.Output_fileS);
    else
        filenames = cellstr('Matrix');
    end
else
    for i = 1:size(matrixAll,3)
        if isfield(cfg,'Output_fileS')
            filenames{1,i} = [cfg.Output_fileS '_' num2str(i)]; %#ok<AGROW>
        else
            filenames{1,i} = ['Matrix_' num2str(i)]; %#ok<AGROW>
        end
    end
end

if isfield(cfg,'Output_filepath') & ~isempty(cfg.Output_filepath)
    [~,~,~,cfg.Output_fileS] = lab_filename(cfg.Output_file);
    if ~isfield(settings.GRAPH,'folder')
        settings.GRAPH.folder = 'GraphResult';
    end
    Output_filepath = fullfile(cfg.Output_filepath,settings.GRAPH.folder);
    warning off %#ok<WNOFF>
    mkdir (Output_filepath);
    warning on %#ok<WNON>
    
    if isfield(settings,'GRAPH') & isfield(settings.GRAPH,'deleteold') & settings.GRAPH.deleteold == true
        warning off %#ok<WNOFF>
        if ~isfield(cfg,'listold')
            if exist(fullfile(Output_filepath,[cfg.Output_fileS '_GRAPH.mat']),'file')
                disp('     delete GraphAnalysis from previous run')
                delete(fullfile(Output_filepath,[cfg.Output_fileS '_GRAPH.mat']));
            end
            cfg.listold = cellstr(fullfile(Output_filepath,[cfg.Output_fileS '_GRAPH.mat']));
        elseif min(~strcmp(cfg.listold,fullfile(Output_filepath,[cfg.Output_fileS '_GRAPH.mat'])))
            if exist(fullfile(Output_filepath,[cfg.Output_fileS '_GRAPH.mat']),'file')
                disp('     delete GraphAnalysis from previous run')
                delete(fullfile(Output_filepath,[cfg.Output_fileS '_GRAPH.mat']));
            end
            cfg.listold = [cfg.listold cellstr(fullfile(Output_filepath,[cfg.Output_fileS '_GRAPH.mat']))];
        end
        warning on %#ok<WNON>
    end
    try %#ok<TRYNC>
        load(fullfile(Output_filepath,[cfg.Output_fileS '_GRAPH.mat']));
    end
    if exist('Result','var') & isfield(Result,'MatrixAll') & ~isempty(Result.MatrixAll)
        if size(matrixAll,1) == size(Result.MatrixAll,1) & ...
                size(matrixAll,2) == size(Result.MatrixAll,2)
            if doIncrease == false
                matrixAll = cat(3,Result.MatrixAll,matrixAll);
                Result.MatrixAll = matrixAll;
            else
                Result.MatrixAll = matrixAll;
            end
        else
            disp ('    Skip matrix, size not consistent to already processed matrices')
            return
        end
        if isfield(Result,'data')
            resultnr = size(Result.data,2) + 1;
        elseif isfield(Result,'datamean')
            resultnr = size(Result.datamean,2) + 1;
        else
            resultnr = 1;
        end
    else
        Result.MatrixAll = matrixAll;
        resultnr = 1;
    end
else
    Result.MatrixAll = matrixAll;
    resultnr = 1;
end

maxmatrixnr = size(matrixAll,3);
if isfield(settings.GRAPH,'maxmatrix') & ~isempty(settings.GRAPH.maxmatrix) & size(matrixAll,3) >= settings.GRAPH.maxmatrix
    % last matrix to calculate, includes calculation of average matrix
    maxmatrixnr = settings.GRAPH.maxmatrix;
end
if isfield(settings.GRAPH,'avgmatrix') & settings.GRAPH.avgmatrix == true & isfield(settings.GRAPH,'numavg')
    if size(matrixAll,3) >= settings.GRAPH.numavg
        avgmatrix = 1;
    elseif isempty(settings.GRAPH.numavg) & (~isfield(cfg,'lastmatrix') | cfg.lastmatrix == true) 
        avgmatrix = 1;
    else
        avgmatrix = 0;
    end
elseif ~isfield(settings.GRAPH,'avgmatrix') & size(matrixAll,3) > 1
    avgmatrix = 1;
else
    avgmatrix = 0;
end

if min(matrixAll(:)) < 0
    disp('    taking absolute values for calculation of MST, eccentricity, betweenneess and efficiency')
end

if ~exist('Result','var') | ~isfield(Result,'results')
    Result.results = {};
    Result.resultsmean = {};
    Result.results2 = {};
end

if isfield(Result,'filenames') & doIncrease == false
    Result.filenames = [Result.filenames filenames];
else
    Result.filenames = filenames;
end

matrixstart = resultnr;
if matrixstart > (maxmatrixnr+avgmatrix)
    return
end
for matrixnr = matrixstart:(maxmatrixnr+avgmatrix)
    if matrixnr == (maxmatrixnr+1)
        disp('      Matrix average');
        if isfield(settings.GRAPH,'numavg') & ~isempty(settings.GRAPH.numavg)
            matrix = mean(matrixAll(:,:,1:settings.GRAPH.numavg),3);
        else
            matrix = mean(matrixAll,3);
        end
    else
        disp(['      Matrix ' num2str(matrixnr)]);
        matrix = matrixAll(:,:,matrixnr);
        Result.ResultNames{1,resultnr} = Result.filenames{1,matrixnr};
    end
    
    % Reduce matrix to Mappings if necessary
    if isfield(settings.GRAPH,'Mappings') & ~isempty(settings.GRAPH.Mappings) & ...
            settings.GRAPH.Mappings.mappingsChannels == size(matrix,1)
        settings.MATRIX.Mappings = settings.GRAPH.Mappings;
        settings.MATRIX.MappingsMode = settings.GRAPH.MappingsMode;
        R = lab_reduce_matrix2mappings(matrix,settings,true);
        matrix = R.matrix;
        clearvars settings R
    end
    numchannels = size(matrix,1);
    matrix(1:numchannels+1:end) = 0;
    Result.matrix = matrix;
    
    if isfield(settings.GRAPH,'exclude') & ~isempty(settings.GRAPH.exclude)
        Idx = setdiff(1:size(matrix,1),settings.GRAPH.exclude);
        matrix = matrix(Idx,Idx);
        numchannels = size(matrix,1);
        matrix(1:numchannels+1:end) = 0;
        Result.matrix = matrix;
    end
    
    % connectivity weights
    matrixtmp = lab_rm_diagonal(matrix);
    Sindex = find(strcmp(Result.results,'Weight'),1);
    Mindex = find(strcmp(Result.resultsmean,'Weight'),1);
    if isempty(Sindex)
        Result.results{end+1,1} = 'Weight';
        Sindex = size(Result.results,1);
        Result.resultsmean{end+1,1} = 'Weight';
        Mindex = size(Result.resultsmean,1);
    end
    Result.data(:,resultnr,Sindex) = mean(matrixtmp,1);
    Result.datamean(Mindex,resultnr) = mean(Result.data(:,resultnr,Sindex));
    clearvars matrixtmp
    
    % rank matrix if necessary
    if isfield(settings,'GRAPH') & isfield(settings.GRAPH,'rankmatrix') & settings.GRAPH.rankmatrix == true
        matrix = lab_rank_matrix(matrix,settings.GRAPH.rankorder);
        matrixtmp = lab_rm_diagonal(matrix);
        Sindex = find(strcmp(Result.results,'WeightRanked'),1);
        Mindex = find(strcmp(Result.resultsmean,'WeightRanked'),1);
        if isempty(Sindex)
            Result.results{end+1,1} = 'WeightRanked';
            Sindex = size(Result.results,1);
            Result.resultsmean{end+1,1} = 'WeightRanked';
            Mindex = size(Result.resultsmean,1);
        end
        Result.data(:,resultnr,Sindex) = mean(matrixtmp,1);
        Result.datamean(Mindex,resultnr) = mean(Result.data(:,resultnr,Sindex));
        clearvars matrixtmp
        isranked = true;
    else
        isranked = false;
    end
    
    % get adjacency matrix
    adjmatrix = zeros(numchannels,numchannels);
    adjmatrix(abs(matrix) > 0) = 1;
    
    % calculate matrixline
    % [tmp1,tmp2,tmp3] = find(matrix);
    % matrixline = [tmp1 tmp2 tmp3];
    % clearvars tmp1 tmp2 tmp3
    
    % calculate inverse matrix
    [matrixI,matrixlineI] = lab_invmatrix(matrix);
    
    fprintf('      do ')
    
    % Clustering Coefficient = average "intensity" of triangles around a node
    if max(strcmp(settings.GRAPH.measure,'ClusteringCoeff'))
        fprintf('ClustCoeff ')
        Sindex = find(strcmp(Result.results,'ClusteringCoeff'),1);
        Mindex = find(strcmp(Result.resultsmean,'AvgClusteringCoeff'),1);
        if isempty(Sindex)
            Result.results{end+1,1} = 'ClusteringCoeff';
            Sindex = size(Result.results,1);
            Result.resultsmean{end+1,1} = 'AvgClusteringCoeff';
            Mindex = size(Result.resultsmean,1);
        end
        if min(matrix(:)) < 0
            matrixtmp = tril(matrix,-1) + tril(matrix,-1)';
            matrixtmp = abs(matrixtmp);
            Result.data(:,resultnr,Sindex)= lab_clusteringcoeff(matrixtmp);
            clearvars matrixtmp
        else
            Result.data(:,resultnr,Sindex)= lab_clusteringcoeff(matrix);
        end
        Result.datamean(Mindex,resultnr) = mean(Result.data(:,resultnr,Sindex));
    end
    
    % Norm of Clustering Coeff
    if max(strcmp(settings.GRAPH.measure,'ClusteringCoeff_Norm')) | max(strcmp(settings.GRAPH.measure,'ShortestPath_Norm'))
        if isranked == true
            NormFile = mfilename('fullpath');
            tmp = strfind(NormFile,filesep);
            if ~isempty(tmp)
                NormFile = NormFile(1:tmp(end));
            else
                NormFile = pwd;
            end
            clearvars tmp
            warning off %#ok<WNOFF>
            mkdir(fullfile(NormFile,'GraphNorm'));
            warning on %#ok<WNON>
            NormFile = fullfile(NormFile,'GraphNorm');
            NormFile = fullfile(NormFile,['GraphNorm_' num2str(size(matrix,1)) '_' ...
                num2str(settings.GRAPH.rankorder) '_' num2str(settings.GRAPH.randnumber) '.mat']);
            if exist(NormFile,'file')
                load(NormFile);
            end
        end
        if max(strcmp(settings.GRAPH.measure,'ClusteringCoeff_Norm'))
            docw = 1;
        else
            docw = 0;
        end
        if max(strcmp(settings.GRAPH.measure,'ShortestPath_Norm'))
            dolambda = 1;
        else
            dolambda = 0;
        end
        if ~exist('cw_r','var') | ~exist('lambda_r','var')
            [cw_r,lambda_r] = lab_graph_randmatrix(matrix,settings.GRAPH.randnumber,settings.GRAPH.randiter,docw,dolambda);
            if isranked == true & exist('NormFile','var')
                save(NormFile,'cw_r','lambda_r');
            end
        end
    end
    if max(strcmp(settings.GRAPH.measure,'ClusteringCoeff_Norm'))
        Mindex = find(strcmp(Result.resultsmean,'AvgClusteringCoeff'),1);
        tmp = Result.datamean(Mindex,resultnr) / cw_r;
        Mindex = find(strcmp(Result.resultsmean,'AvgClusteringCoeff_Norm'), 1);
        if isempty(Mindex)
            Result.resultsmean{end+1,1} = 'AvgClusteringCoeff_Norm';
            Mindex = size(Result.resultsmean,1);
        end
        Result.datamean(Mindex,resultnr) = tmp;
        clearvars cw_r tmp
    end
    
    % Shortest Paths
    if max(strcmp(settings.GRAPH.measure,'ShortestPath'))
        fprintf('ShortPath ')
        Sindex = find(strcmp(Result.results,'ShortestPath'),1);
        Mindex = find(strcmp(Result.resultsmean,'AvgShortestPath'),1);
        Tindex = find(strcmp(Result.results2,'ShortestPath'),1);
        if isempty(Sindex)
            Result.results{end+1,1} = 'ShortestPath';
            Sindex = size(Result.results,1);
            Result.resultsmean{end+1,1} = 'AvgShortestPath';
            Mindex = size(Result.resultsmean,1);
            Result.results2{end+1,1} = 'ShortestPath';
            Tindex = size(Result.results2,1);
            Result.results2{end+1,1} = 'ShortestPathDistances';
            Result.results2{end+1,1} = 'ShortestPathNumber';
        end
        [distance, numberpaths]=distance_wei(matrixI);
        Result.data2{Tindex+1,resultnr} = distance;
        Result.data2{Tindex+2,resultnr} = numberpaths;
        distance(1:size(distance,1)+1:end) = Inf;
        distance = distance.^-1;
        Result.data(:,resultnr,Sindex) = (sum(distance,2) / (size(distance,1)-1)).^-1;
        Result.datamean(Mindex,resultnr) = (sum(distance(:)) / (size(distance,1)*(size(distance,1)-1))).^-1;
        if max(strcmp(settings.GRAPH.measure,'ShortestPath_Norm'))
            tmp = Result.datamean(Mindex,resultnr) / lambda_r;
            Mindex = find(strcmp(Result.resultsmean,'ShortestPath_Norm'), 1);
            if isempty(Mindex)
                Result.resultsmean{end+1,1} = 'ShortestPath_Norm';
                Mindex = size(Result.resultsmean,1);
            end
            Result.datamean(Mindex,resultnr) = tmp;
            clearvars tmp lambda_r
        end
        clearvars distance numberpaths 
        [~,Result.data2{Tindex,resultnr}] = dijkstra(adjmatrix,matrixI);
    end
    
    % Degree = number of links connected to the node
    if max(strcmp(settings.GRAPH.measure,'Degree'))
        fprintf('Degree ')
        Sindex = find(strcmp(Result.results,'Degree'), 1);
        if isempty(Sindex)
            Result.results{end+1,1} = 'Degree';
            Sindex = size(Result.results,1);
        end
        if size(unique(matrix),1) == 2
            Result.data(:,resultnr,Sindex)  = degrees_und(matrix);
        else
            matrixtmp = lab_rm_diagonal(matrix);
            Result.data(:,resultnr,Sindex)  = mean(matrixtmp,1);
        end
    end
    
    % Degree correlation
    if max(strcmp(settings.GRAPH.measure,'DegreeCorrelation'))
        fprintf('DegreeCorrelation ')
        Mindex = find(strcmp(Result.resultsmean,'DegreeCorrelation'),1);
        if isempty(Mindex)
            Result.resultsmean{end+1,1} = 'DegreeCorrelation';
            Mindex = size(Result.resultsmean,1);
        end
        matrixtmp = matrix;
        matrixtmp(1:size(matrixtmp,1)+1:end) = 0;
        degree = sum(matrixtmp,1);
        matrixtmp = (tril(matrixtmp) + triu(matrixtmp)') ./ 2;
        [tmp1,tmp2] = find(matrixtmp>0);
        Result.datamean(Mindex,resultnr) = corr(degree(tmp1)',degree(tmp2)');
        clearvars matrixtmp degree
    end
    
    % Degree diversity
    if max(strcmp(settings.GRAPH.measure,'DegreeDiversity'))
        fprintf('DegreeDiversity ')
        Mindex = find(strcmp(Result.resultsmean,'DegreeDiversity'),1);
        if isempty(Mindex)
            Result.resultsmean{end+1,1} = 'DegreeDiversity';
            Mindex = size(Result.resultsmean,1);
        end
        matrixtmp = matrix;
        matrixtmp(1:size(matrixtmp,1)+1:end) = 0;
        degree = sum(matrixtmp,1);
        Result.datamean(Mindex,resultnr) = mean(degree.^2)/mean(degree);
        clearvars matrixtmp degree
    end
    
    % Density = fraction of present connections to possible connections
    if max(strcmp(settings.GRAPH.measure,'Density'))
        fprintf('Density ')
        Mindex = find(strcmp(Result.resultsmean,'Density'), 1);
        if isempty(Mindex)
            Result.resultsmean{end+1,1} = 'Density';
            Mindex = size(Result.resultsmean,1);
        end
        if min(matrix(:)) < 0
            Result.datamean(Mindex,resultnr) = density_und(matrix);
        else
            Result.datamean(Mindex,resultnr) = density_dir(matrix);
        end
    end
    
    % Eccentricity
    if max(strcmp(settings.GRAPH.measure,'Eccentricity'))
        fprintf('Eccentricity ')
        Sindex = find(strcmp(Result.results,'Eccentricity'),1);
        Mindex = find(strcmp(Result.resultsmean,'MaxEccentricity'),1);
        Tindex = find(strcmp(Result.results2,'EccentricityCentralNodes'),1);
        if isempty(Sindex)
            Result.results{end+1,1} = 'Eccentricity';
            Sindex = size(Result.results,1);
            Result.resultsmean{end+1,1} = 'MaxEccentricity';
            Mindex = size(Result.resultsmean,1);
            Result.resultsmean{end+1,1} = 'Radius';
            Result.resultsmean{end+1,1} = 'Diameter';
            Result.results2{end+1,1} = 'EccentricityCentralNodes';
            Tindex = size(Result.results2,1);
            Result.results2{end+1,1} = 'EccentricityPeripheralNodes';
        end
        [Result.data(:,resultnr,Sindex),Result.datamean(Mindex+1,resultnr), ...
            Result.datamean(Mindex+2,resultnr),CenterNodes,PeripheralNodes] = grEccentricity(matrixlineI);
        Result.data2{Tindex,resultnr} = num2str(CenterNodes);
        Result.data2{Tindex+1,resultnr}  = num2str(PeripheralNodes);
        Result.datamean(Mindex,resultnr) = max(Result.data(:,resultnr,Sindex));
        clearvars CenterNodes PeripheralNodes
    end
    
    % Edge/Node betweenness centrality = fraction of all shortest paths
    % in the network that contain a given edge
    if max(strcmp(settings.GRAPH.measure,'Betweenness'))
        fprintf('Betweenness ')
        Sindex = find(strcmp(Result.results,'Betweenness'),1);
        Mindex = find(strcmp(Result.resultsmean,'MaxBetweenness'),1);
        Tindex = find(strcmp(Result.results2,'EdgeBetweenness'),1);
        if isempty(Sindex)
            Result.results{end+1,1} = 'Betweenness';
            Sindex = size(Result.results,1);
            Result.resultsmean{end+1,1} = 'MaxBetweenness';
            Mindex = size(Result.resultsmean,1);
            Result.results2{end+1,1} = 'EdgeBetweenness';
            Tindex = size(Result.results2,1);
        end
        [Result.data2{Tindex,resultnr},NodeBetweenness] = edge_betweenness_wei(matrixI);
        Result.data(:,resultnr,Sindex) = NodeBetweenness' / ((numchannels-1)*(numchannels-2));
        Result.datamean(Mindex,resultnr) = max(Result.data(:,resultnr,Sindex));
        clearvars NodeBetweenness
    end
    
    % Efficiency = average of inverse shortest path length
    if max(strcmp(settings.GRAPH.measure,'Efficiency'))
        fprintf('Efficiency ')
        Mindex = find(strcmp(Result.resultsmean,'Efficiency'), 1);
        if isempty(Mindex)
            Result.resultsmean{end+1,1} = 'Efficiency';
            Mindex = size(Result.resultsmean,1);
        end 
        Result.datamean(Mindex,resultnr) = efficiency_wei(abs(matrix));
    end
    
    % Eigenvector centrality
    if max(strcmp(settings.GRAPH.measure,'EigenvectorCentrality'))
        fprintf('Eigenvector ')
        Sindex = find(strcmp(Result.results,'EigenvectorCentrality'),1);
        Mindex = find(strcmp(Result.resultsmean,'MaximalEigenvalue'),1);
        if isempty(Sindex)
            Result.results{end+1,1} = 'EigenvectorCentrality';
            Sindex = size(Result.results,1);
            Result.resultsmean{end+1,1} = 'MaximalEigenvalue';
            Mindex = size(Result.resultsmean,1);
            Result.resultsmean{end+1,1} = 'NodeMaxEigenvectorCentrality';
        end
        [EC,EM] = eigenvector_centrality_und(matrix);
        Result.data(:,resultnr,Sindex) = EC';
        Result.datamean(Mindex,resultnr) = EM;
        Result.datamean(Mindex+1,resultnr) = find(EC == max(EC),1);
        clearvars EC EM
    end
    
    % Transitivity
    if max(strcmp(settings.GRAPH.measure,'Transitivity'))
        fprintf('Transitivity ')
        Mindex = find(strcmp(Result.resultsmean,'Transitivity'),1);
        if isempty(Mindex)
            Result.resultsmean{end+1,1} = 'Transitivity';
            Mindex = size(Result.resultsmean,1);
        end
        if min(matrix(:)) >= 0
            Result.datamean(Mindex,resultnr) = transitivity_wu(matrix);
        else
            Result.datamean(Mindex,resultnr) = transitivity_wd(matrix);
        end
    end
    
    % Modularity = statistic that quantifies the degree to which the
    % network may be subdivided into such clearly delineated groups
    if max(strcmp(settings.GRAPH.measure,'Modularity'))
        fprintf('Modularity ')
        Sindex = find(strcmp(Result.results,'ModularityCS'),1);
        Mindex = find(strcmp(Result.resultsmean,'Modularity'),1);
        if isempty(Sindex)
            Result.results{end+1,1} = 'ModularityCS';
            Sindex = size(Result.results,1);
            Result.resultsmean{end+1,1} = 'Modularity';
            Mindex = size(Result.resultsmean,1);
        end
        if min(matrix(:)) < 0
            [ModularityCS, Result.datamean(Mindex,resultnr)]=modularity_und(matrix);
        else
            [ModularityCS, Result.datamean(Mindex,resultnr)]=modularity_dir(matrix);
        end
        Result.data(:,resultnr,Sindex) = ModularityCS';
        clearvars ModularityCS
    end
    
    % Matrix Ref
    if max(strcmp(settings.GRAPH.measure,'Matrix-Ref')) & isfield(settings.GRAPH,'MatrixRef') & ~isempty(settings.GRAPH.MatrixRef) & ...
            size(matrix,1) == size(settings.GRAPH.MatrixRef,1) & size(matrix,2) == size(settings.GRAPH.MatrixRef,2)
        fprintf('Matrix-Ref ')
        Mindex = find(strcmp(Result.resultsmean,'Matrix-Ref'),1);
        if isempty(Mindex)
            Result.resultsmean{end+1,1} = 'Matrix-Ref';
            Mindex = size(Result.resultsmean,1);
        end
        Result.datamean(Mindex,resultnr) = corr(matrix(:),settings.GRAPH.MatrixRef(:),'type','Kendall');
        % Result.datamean(Mindex,resultnr) = corr(matrix(:),settings.GRAPH.MatrixRef(:));
    end
    
    if max(~cellfun(@isempty,strfind(settings.GRAPH.measure,'MST-')))
        fprintf('MST ')
        %------------------------------
        % MST
        %------------------------------
        [MSTmatrix,Result.MSTweight(1,resultnr),MSTline] = lab_MST(matrix);
        Result.MST(:,:,resultnr) = MSTmatrix;
        
        % Degree
        if max(strcmp(settings.GRAPH.measure,'MST-Degree'))
            fprintf('MST-Degree ')
            Sindex = find(strcmp(Result.results,'MST-Degree'),1);
            if isempty(Sindex)
                Result.results{end+1,1} = 'MST-Degree';
                Sindex = size(Result.results,1);
            end
            [~,~,Result.data(:,resultnr,Sindex)] = degrees_dir(MSTmatrix);
        end
        
        % Degree correlation
        if max(strcmp(settings.GRAPH.measure,'MST-DegreeCorrelation'))
            fprintf('MST-DegreeCorrelation ')
            Mindex = find(strcmp(Result.resultsmean,'MST-DegreeCorrelation'),1);
            if isempty(Mindex)
                Result.resultsmean{end+1,1} = 'MST-DegreeCorrelation';
                Mindex = size(Result.resultsmean,1);
            end
            degree = sum(MSTmatrix,1);
            matrixtmp = tril(MSTmatrix);
            [tmp1,tmp2] = find(matrixtmp>0);
            Result.datamean(Mindex,resultnr) = corr(degree(tmp1)',degree(tmp2)');
            clearvars matrixtmp degree
        end
        
        % Degree diversity
        if max(strcmp(settings.GRAPH.measure,'MST-DegreeDiversity'))
            fprintf('MST-DegreeDiversity ')
            Mindex = find(strcmp(Result.resultsmean,'MST-DegreeDiversity'),1);
            if isempty(Mindex)
                Result.resultsmean{end+1,1} = 'MST-DegreeDiversity';
                Mindex = size(Result.resultsmean,1);
            end
            degree = sum(MSTmatrix,1);
            Result.datamean(Mindex,resultnr) = mean(degree.^2)/mean(degree);
            clearvars degree
        end
        
        % Shortest paths
        if max(strcmp(settings.GRAPH.measure,'MST-ShortestPath'))
            fprintf('MST-ShortestPath ')
            Mindex = find(strcmp(Result.resultsmean,'MST-AvgShortestPath'),1);
            Tindex = find(strcmp(Result.results2,'MST-ShortestPaths'),1);
            if isempty(Mindex)
                Result.resultsmean{end+1,1} = 'MST-AvgShortestPath';
                Mindex = size(Result.resultsmean,1);
                Result.results2{end+1,1} = 'MST-ShortestPaths';
                Tindex = size(Result.results2,1);
            end
            distance = distance_bin(MSTmatrix);
            Result.data2{Tindex,resultnr} = distance;
            Result.datamean(Mindex,resultnr) = sum(sum(distance)) / (numchannels * (numchannels - 1));
            clearvars distance
        end
        
        % Eccentricity
        if max(strcmp(settings.GRAPH.measure,'MST-Eccentricity'))
            fprintf('MST-Eccentricity ')
            Sindex = find(strcmp(Result.results,'MST-Eccentricity'),1);
            Mindex = find(strcmp(Result.resultsmean,'MST-MaxEccentricity'),1);
            Tindex = find(strcmp(Result.results2,'MST-EccentricityCentralNodes'),1);
            if isempty(Sindex)
                Result.results{end+1,1} = 'MST-Eccentricity';
                Sindex = size(Result.results,1);
                Result.resultsmean{end+1,1} = 'MST-MaxEccentricity';
                Mindex = size(Result.resultsmean,1);
                Result.resultsmean{end+1,1} = 'MST-Radius';
                Result.resultsmean{end+1,1} = 'MST-Diameter';
                Result.results2{end+1,1} = 'MST-EccentricityCentralNodes';
                Tindex = size(Result.results2,1);
                Result.results2{end+1,1} = 'MST-EccentricityPeripheralNodes';
            end
            [Result.data(:,resultnr,Sindex),Result.datamean(Mindex+1,resultnr), ...
                Result.datamean(Mindex+2,resultnr), ...
                CenterNodes,PeripheralNodes] = grEccentricity(MSTline);
            Result.data2{Tindex,resultnr} = num2str(CenterNodes);
            Result.data2{Tindex+1,resultnr} = num2str(PeripheralNodes);
            Result.datamean(Mindex,resultnr) = max(Result.data(:,resultnr,Sindex));
            clearvars CenterNodes PeripheralNodes
        end
        
        % Betweenness
        if max(strcmp(settings.GRAPH.measure,'MST-Betweenness'))
            fprintf('MST-Betweenness ')
            Sindex = find(strcmp(Result.results,'MST-Betweenness'),1);
            Mindex = find(strcmp(Result.resultsmean,'MST-MaxBetweenness'),1);
            Tindex = find(strcmp(Result.results2,'MST-EdgeBetweenness'),1);
            if isempty(Sindex)
                Result.results{end+1,1} = 'MST-Betweenness';
                Sindex = size(Result.results,1);
                Result.resultsmean{end+1,1} = 'MST-MaxBetweenness';
                Mindex = size(Result.resultsmean,1);
                Result.results2{end+1,1} = 'MST-EdgeBetweenness';
                Tindex = size(Result.results2,1);
            end
            [Result.data2{Tindex,resultnr},NodeBetweenness] = edge_betweenness_bin(MSTmatrix);
            Result.data(:,resultnr,Sindex) = NodeBetweenness' / ((numchannels-1)*(numchannels-2));
            Result.datamean(Mindex,resultnr) = max(Result.data(:,resultnr,Sindex));
            clearvars NodeBetweenness
        end
        
        % Efficiency
        if max(strcmp(settings.GRAPH.measure,'MST-Efficiency'))
            fprintf('MST-Efficiency ')
            Mindex = find(strcmp(Result.resultsmean,'MST-Efficiency'),1);
            if isempty(Mindex)
                Result.resultsmean{end+1,1} = 'MST-Efficiency';
                Mindex = size(Result.resultsmean,1);
            end
            Result.datamean(Mindex,resultnr) = efficiency_bin(MSTmatrix);
        end
        
        % Eigenvector centrality
        if max(strcmp(settings.GRAPH.measure,'MST-EigenvectorCentrality'))
            fprintf('MST-EC ')
            Sindex = find(strcmp(Result.results,'MST-EigenvectorCentrality'),1);
            Mindex = find(strcmp(Result.resultsmean,'MST-MaximalEigenvalue'),1);
            if isempty(Sindex)
                Result.results{end+1,1} = 'MST-EigenvectorCentrality';
                Sindex = size(Result.results,1);
                Result.resultsmean{end+1,1} = 'MST-MaximalEigenvalue';
                Mindex = size(Result.resultsmean,1);
                Result.resultsmean{end+1,1} = 'MST-NodeMaxEigenvectorCentrality';
            end
            [EC,EM] = eigenvector_centrality_und(MSTmatrix);
            Result.data(:,resultnr,Sindex) = EC';
            Result.datamean(Mindex,resultnr) = EM;
            Result.datamean(Mindex+1,resultnr) = find(EC == max(EC));
            clearvars EM EC
        end
        
        % MST leave nodes
        if max(strcmp(settings.GRAPH.measure,'MST-LeaveNodes'))
            fprintf('MST-LeaveNodes ')
            Mindex = find(strcmp(Result.resultsmean,'MST-LeaveNodes'),1);
            if isempty(Mindex)
                Result.resultsmean{end+1,1} = 'MST-LeaveNodes';
                Mindex = size(Result.resultsmean,1);
            end
            tmp = sum(MSTmatrix,1);
            Result.datamean(Mindex,resultnr) = sum(tmp == 1);
        end
        
        % MST Ref
        if max(strcmp(settings.GRAPH.measure,'MST-Ref')) & isfield(settings.GRAPH,'MSTref') & ~isempty(settings.GRAPH.MSTref) & ...
                size(MSTmatrix,1) == size(settings.GRAPH.MSTref,1) & size(MSTmatrix,2) == size(settings.GRAPH.MSTref,2)
            fprintf('MST-Ref ')
            Mindex = find(strcmp(Result.resultsmean,'MST-Ref'),1);
            if isempty(Mindex)
                Result.resultsmean{end+1,1} = 'MST-Ref';
                Mindex = size(Result.resultsmean,1);
            end
            tmp = MSTmatrix(settings.GRAPH.MSTref == 1);
            Result.datamean(Mindex,resultnr) = sum(tmp) / length(tmp);
            clearvars tmp
        end
        
        clearvars MSTline MSTmatrix
    end
    
    clearvars matrix matrixI matrixline matrixlineI adjmatrix
    disp('.')
    
    resultnr = resultnr + 1;
end

if avgmatrix == 1
    % Store AVG results in separate variable
    Result.dataAVG = Result.data(:,end,:);
    Result.datameanAVG = Result.datamean(:,end);
    if isfield(Result,'data2')
        Result.data2AVG = Result.data2(:,end);
    end
    if size(Result.data,2) > 1
        Result.data = Result.data(:,1:end-1,:);
        Result.datamean = Result.datamean(:,1:end-1);
        if isfield(Result,'data2')
            Result.data2 = Result.data2(:,1:end-1);
        end
    else
        Result.data = [];
        Result.datamean = [];
        Result.data2 = [];
    end
else
    Result.matrix = mean(matrixAll,3);
    if isfield(settings.GRAPH,'Mappings') & ~isempty(settings.GRAPH.Mappings) & ...
            settings.GRAPH.Mappings.mappingsChannels == size(matrixAll,1)
        settings.MATRIX.Mappings = settings.GRAPH.Mappings;
        settings.MATRIX.MappingsMode = settings.GRAPH.MappingsMode;
        R = lab_reduce_matrix2mappings(Result.matrix,settings,true);
        Result.matrix = R.matrix;
        clearvars settings R
    end
end
if exist('timestamp','var') & ~isempty(timestamp)
    if isfield(Result,'timestamp') & doIncrease == false
        Result.timestamp = [Result.timestamp timestamp(:)'];
    else
        Result.timestamp = timestamp(:)';
    end
end

if isfield(cfg,'Output_filepath') & ~isempty(cfg.Output_filepath)
    % write mean xls
    disp('      Write results to xls');
    if isfield(Result,'datamean') & ~isempty(Result.datamean)
        headertmp = [cellstr('') Result.ResultNames];
        xlsout = [Result.resultsmean num2cell(Result.datamean)];
        xlsout = cat(1,headertmp,xlsout);
        if size(xlsout,2) > 255
            fileout = fullfile(Output_filepath,[cfg.Output_fileS '_mean.xlsx']);
        else
            fileout = fullfile(Output_filepath,[cfg.Output_fileS '_mean.xls']);
        end
        lab_write_xls(fileout,xlsout);
        clearvars xlsout fileout headertmp
    end
    
    % write single variables
    if isfield(Result,'data') & ~isempty(Result.data)
        headertmp = [cellstr('') Result.ResultNames];
        dataall = [];
        for i = 1:size(Result.results,1)
            for j = 1:numchannels
                if isfield(header,'channels') & size(header.channels,1) >= j
                    resultstmp{j,i} = [Result.results{i,1} '_' regexprep(header.channels(j,:),' ','')]; %#ok<AGROW>
                else
                    resultstmp{j,i} = [Result.results{i,1} '_ch' num2str(j)]; %#ok<AGROW>
                end
            end
            dataall = cat(1,dataall,Result.data(:,:,i));
        end
        resultstmp = resultstmp(:);
        xlsout = [resultstmp num2cell(dataall)];
        xlsout = cat(1,headertmp,xlsout);
        xlsout{1,1} = ['C' num2str(size(Result.data,1)) ' R0'];
        if size(xlsout,2) > 255
            fileout = fullfile(Output_filepath,[cfg.Output_fileS '.xlsx']);
        else
            fileout = fullfile(Output_filepath,[cfg.Output_fileS '.xls']);
        end
        lab_write_xls(fileout,xlsout);
        clearvars xlsout dataall fileout resultstmp headertmp i j
    end
    
    if isfield(Result,'datameanAVG')
        % write AVG results
        headertmp = {'','AVG'};
        xlsout = [Result.resultsmean num2cell(Result.datameanAVG)];
        xlsout = cat(1,headertmp,xlsout);
        if size(xlsout,2) > 255
            fileout = fullfile(Output_filepath,[cfg.Output_fileS '_AVGmean.xlsx']);
        else
            fileout = fullfile(Output_filepath,[cfg.Output_fileS '_AVGmean.xls']);
        end
        lab_write_xls(fileout,xlsout);
        clearvars xlsout fileout headertmp
    end
        
    % write AVG single
    if isfield(Result,'dataAVG')
        headertmp = {'','AVG'};
        dataallAVG = [];
        resultstmp = {};
        for i = 1:size(Result.results,1)
            for j = 1:numchannels
                if isfield(header,'channels') & size(header.channels,1) >= j
                    resultstmp{j,i} = [Result.results{i,1} '_' regexprep(header.channels(j,:),' ','')]; %#ok<AGROW>
                else
                    resultstmp{j,i} = [Result.results{i,1} '_ch' num2str(j)]; %#ok<AGROW>
                end
            end
            dataallAVG = cat(1,dataallAVG,Result.dataAVG(:,:,i));
        end
        resultstmp = resultstmp(:);         
        xlsout = [resultstmp num2cell(dataallAVG)];
        xlsout = cat(1,headertmp,xlsout);
        xlsout{1,1} = ['C' num2str(size(Result.dataAVG,1)) ' R0'];
        if size(xlsout,2) > 255
            fileout = fullfile(Output_filepath,[cfg.Output_fileS '_AVG.xlsx']);
        else
            fileout = fullfile(Output_filepath,[cfg.Output_fileS '_AVG.xls']);
        end
        lab_write_xls(fileout,xlsout);
        clearvars headertmp resultstmp dataallAVG xlsout fileout i j
    end
    
    
    % Write Images and matrixes
    f = figure('Visible','off');
    colormap('gray');
    cmap = colormap;
    colormap(flipud(cmap));
    imagesc(Result.matrix);
    lab_print_figure(fullfile(Output_filepath,[cfg.Output_fileS '.tif']),f);
    clf;
    if isfield(Result,'MST')
        imagesc(Result.MST(:,:,end));
        lab_print_figure(fullfile(Output_filepath,[cfg.Output_fileS '_MST.tif']),f);
        clf;
    end
    close;
    % write matrix
    matrixfileout = fullfile(Output_filepath,[cfg.Output_fileS '_matrix.txt']);
    if exist(matrixfileout,'file')
        delete(matrixfileout);
    end
    dlmwrite(matrixfileout,Result.matrix,'delimiter','\t','precision', 6);
    % write MST-matrix
    if isfield(Result,'MST')
        matrixfileout = fullfile(Output_filepath,[cfg.Output_fileS '_MST.txt']);
        if exist(matrixfileout,'file')
            delete(matrixfileout);
        end
        dlmwrite(matrixfileout,Result.MST(:,:,end),'delimiter','\t','precision', 6);
    end
    
    if isfield(header,'patient')
        Result.patient = header.patient;
    elseif isfield(cfg,'patient')
        Result.patient = cfg.patient;
    end
    if isfield(header,'freqband') & ~isempty(header.freqband)
        Result.freqband = ['F' num2str(header.freqband(1)) 'F' num2str(header.freqband(2))];
    end
    save(fullfile(Output_filepath,[cfg.Output_fileS '_GRAPH.mat']),'Result');
end

return
