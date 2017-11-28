% Calculate kmeans clustering for matrices
%
% Result = lab_calculate_kmeans(matrix,maxclusters)
%
%    matrix        array (nchans x nchans x trials)
%    maxclusters   maximal number of clusters for kmeans-clustering
%
% written by F. Hatz 2012

function Result = lab_calculate_kmeans(matrix,maxclusters,Output_filepath,WriteMatrices,doPCA,nofigs)
disp('   Calculate Kmeans-clustering of matrixes')

if ~exist('nofigs','var')
    nofigs = false;
end

if ~exist('matrix','var')
    [settings,skipprocessing] = lab_set_calculate_kmeans;
    if skipprocessing == 1 | ~exist(settings.searchfolder,'dir')
        Result = [];
        return
    end
    maxclusters = settings.maxclusters;
        
    cfg.SEARCH = settings;
    calc = lab_search_files(cfg);
    Filelist = calc.Filelist;
    clearvars cfg calc
    matrix = [];
    for i = 1:size(Filelist,2)
        matrixtmp = lab_read_data(Filelist{1,i});
        if size(matrixtmp,1) == size(matrixtmp,2)
            if isempty(matrix)
                matrix = matrixtmp;
            elseif size(matrixtmp,1) == size(matrix,1) & size(matrixtmp,2) == size(matrix,2)
                matrix = cat(3,matrix,matrixtmp);
            end
        end
    end
    if isempty(matrix)
        Result = [];
        return
    end
elseif ~exist('maxclusters','var')
    [settings,skipprocessing] = lab_set_calculate_kmeans([],1);
    if skipprocessing == 1
        Result = [];
        return
    end
    maxclusters = settings.maxclusters;
end

if ~exist('maxclusters','var') | isempty(maxclusters)
    maxclusters = 120;
end

if maxclusters > size(matrix,3)
    maxclusters = size(matrix,3);
end

if ~exist('WriteMatrices','var')
    if exist('settings','var') & isfield(settings,'WriteMatrices')
        WriteMatrices = settings.WriteMatrices;
    else
        WriteMatrices = true;
    end
end
if ~exist('doPCA','var')
    if exist('settings','var') & isfield(settings,'doPCA')
        doPCA = settings.doPCA;
    else
        doPCA = false;
    end
end
if ~exist('Output_filepath','var') & exist('settings','var') & isfield(settings,'searchfolder')
    Output_filepath = fullfile(settings.searchfolder,'KmeansClustering');
    warning off %#ok<WNOFF>
    mkdir(Output_filepath)
    warning on %#ok<WNON>
end

% K-means
if size(matrix,3) > 20
    if doPCA == true
        ResultPCA = lab_matrix2clustering(matrix,1);
        matrixV = ResultPCA.PCA.score(:,1:ResultPCA.nFactors);
    else
        matrixV = reshape(matrix,size(matrix,1)*size(matrix,2),size(matrix,3))';
    end
    for i = 2:maxclusters
        disp(['     Calculate cluster ' num2str(i) ' of ' num2str(maxclusters) ' clusters'])
        [Result.clusters(:,i),Result.centroid{i},Result.sumclusters{i}] = kmeans(matrixV,i);
    end
    for i = 1:maxclusters-1
        for j = 1:i+1
            selection = (Result.clusters(:,i+1) == j);
            if sum(selection) > 1
                tmp(j,:) = mean(matrixV(Result.clusters(:,i+1) == j,:));
            else
                tmp(j,:) = matrixV(Result.clusters(:,i+1) == j,:);
            end
        end
        Result.meanvar(1,i) = mean(var(tmp));
        Result.meanvarC(1,i) = mean(var(tmp)) /  mean(var(matrixV));
        clear tmp
    end
    if exist('Output_filepath','var') & exist(Output_filepath,'dir')
        % Plot normalized variance
        save(fullfile(Output_filepath,'KmeansClustering.mat'),'Result');
        if nofigs == true
            fig1 = figure('Color',[1 1 1],'NumberTitle','off','Name','Kmeans-Clustering -- Normalized variance','Visible','off');
        else
            fig1 = figure('Color',[1 1 1],'NumberTitle','off','Name','Kmeans-Clustering -- Normalized variance','Menubar','none');
            m1 = uimenu(fig1,'Label','File');
            uimenu(m1,'Label','Save Picture','Callback',@(~,~)lab_print_figure);
            uimenu(m1,'Label','Copy Picture','Callback','print -dmeta;');
            uimenu(m1,'Label','Close','Callback','close;');
        end
        plot(Result.meanvarC');
        title('Kmeans-Clustering -- Normalized variance')
        xlabel('Number of clusters');
        lab_print_figure(fullfile(Output_filepath,'KmeansClustering.tif'),fig1);
        if nofigs == true
            close(fig1);
        end
        
        % Plot Cluster-Allocation as matrix
        if nofigs == true
            fig2 = figure('Color',[1 1 1],'NumberTitle','off','Name','Kmeans-Clustering','Visible','off');
        else
            fig2 = figure('Color',[1 1 1],'NumberTitle','off','Name','Kmeans-Clustering','Menubar','none');
            m2 = uimenu(fig2,'Label','File');
            uimenu(m2,'Label','Save Picture','Callback',@(~,~)lab_print_figure);
            uimenu(m2,'Label','Copy Picture','Callback','print -dmeta;');
            uimenu(m2,'Label','Close','Callback','close;');
        end
        clusters = Result.clusters';
        clusters(2:end,:) = clusters(2:end,:) - 1;
        for i = 1:size(clusters,1)
            tmp = max(clusters(:)) / max(clusters(i,:));
            clusters(i,:) = floor(clusters(i,:) * tmp);
            clearvars tmp
        end
        imagesc(clusters);
        clearvars clusters
        title('Kmeans-Clustering')
        xlabel('Matrices');
        ylabel('Number of clusters');
        lab_print_figure(fullfile(Output_filepath,'KmeansClustering_Clusters.tif'),fig2);
        if nofigs == true
            close(fig2);
        end
        
        % Write normalized variance of clustering
        xlsout = {'Number of clusters','Variance','Normalized variance'};
        xlsout = cat(1,xlsout,[num2cell((2:maxclusters)') num2cell(Result.meanvar') num2cell(Result.meanvarC')]);
        lab_write_xls(fullfile(Output_filepath,'KmeansClustering.xlsx'),xlsout);
        
        % Write cluster-allocation of matrices
        lab_write_xls(fullfile(Output_filepath,'KmeansClustering_Clusters.xlsx'),num2cell(Result.clusters));
        
        % Write mean matrices of every cluster
        if WriteMatrices == true
            disp('     write resulting matrices for all number of clusters')
            for i = 1:maxclusters-1
                Matrix_filepath = fullfile(Output_filepath,['Clusters_' num2str(i+1)]);
                warning off %#ok<WNOFF>
                mkdir(Matrix_filepath);
                warning on %#ok<WNON>
                index = Result.clusters(:,i+1);
                for j = 1:i+1
                    fig1 = lab_plot_matrix(matrix(:,:,index == j),1);
                    lab_print_figure(fullfile(Matrix_filepath,['Matrix_' num2str(j) '.jpg']),fig1);
                    close(fig1)
                    lab_write_matrix(fullfile(Matrix_filepath,[num2str(j) '_matrix.txt']),mean(matrix(:,:,index == j),3));
                end
            end
        end
    end
else
    disp('   Abort (need more matrixes over time)')
    Result.meanvarC = [];
    Result.meanvar = [];
end