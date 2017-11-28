function Output = lab_matrix_clustering(Result,settings)

Output = [];

if ~exist('Result','var')
    [tmp1,tmp2] = uigetfile('*.mat','Select Matlab-Container with results of ''Collect Connectivity''');
    if isnumeric(tmp1)
        return
    end
    load(fullfile(tmp2,tmp1));
    if ~exist('Result','var')
        disp('Abort: no valid Matlab-Container')
        return
    end
    clearvars tmp1 tmp2
end

if ~exist('settings','var') | ~isfield(settings,'maxclusters')
    [settings,skipprocessing] = lab_set_calculate_kmeans([],1);
    if skipprocessing == 1
        return
    end
end
doPCA = settings.doPCA;
WriteMatrices = settings.WriteMatrices;
maxclusters = settings.maxclusters;

if ~isstruct(Result)
    disp('Abort: wrong input data')
    return
end
if isfield(Result,'matrix')
    dofreqs = false;
    Freqs = NaN;
else
    Freqs = fieldnames(Result);
    flags = false(1,length(Freqs));
    for i = 1:length(Freqs)
        if isstruct(Result.(Freqs{i})) & isfield(Result.(Freqs{i}),'matrix')
            flags(i) = true;
        end
    end
    Freqs = Freqs(flags);
    dofreqs = true;
    clearvars i flags
end
for F = 1:length(Freqs);
    if dofreqs == true
        Rtmp = Result.(Freqs{F});
        Foldername = Freqs{F};
    else
        Rtmp = Result;
        Foldername = 'Kmeans';
    end
    Vars = fieldnames(Rtmp.matrix);
    flags = false(1,length(Vars));
    for V = 1:length(Vars)
        if isnumeric(Rtmp.matrix.(Vars{V}))& length(Rtmp.matrix.(Vars{V})) > 1 & ...
                size(Rtmp.matrix.(Vars{V}),1) == size(Rtmp.matrix.(Vars{V}),2)
            flags(V) = true;
        end
    end
    Vars = Vars(flags);
    clearvars flags
    for V = 1:length(Vars)
        Output_filepath = fullfile(settings.searchfolder,[Foldername '_' Vars{V}]);
        warning off %#ok<WNOFF>
        mkdir(Output_filepath);
        warning on %#ok<WNON>
        matrixtmp = Rtmp.matrix.(Vars{V});
        if isfield(Rtmp,'subjects') & isfield(Rtmp.subjects,Vars{V})
            subjects = Rtmp.subjects.(Vars{V});
            if ~exist('results','var')
                if ~isfield(settings,'results') | ~exist(settings.results,'file')
                    [tmp1,tmp2] = uigetfile('*.xlsx','Select file with results');
                    settings.results = fullfile(tmp2,tmp1);
                    clearvars tmp1 tmp2
                end
                results = lab_read_xls(settings.results);
            end
            result = lab_get_group(subjects,results);
            if length(result) == size(Rtmp.matrix.(Vars{V}),3)
                index = result>0;
                if ~isempty(index)
                    matrixtmp = matrixtmp(:,:,index);
                    result = result(index);
                    subjects = subjects(index);
                end
            end
        end
        Output.(Foldername).(Vars{V}) = lab_calculate_kmeans(matrixtmp,maxclusters,Output_filepath,WriteMatrices,doPCA,true);
        Clusters = Output.(Foldername).(Vars{V}).clusters;
        if exist('subjects','var')
            Clusters(:,1) = result(:);
            header = cellstr(num2str((2:size(Clusters,2))'))';
            header = cat(2,cellstr('group'),header);
            xlsout = cat(1,header,num2cell(Clusters));
            subjects = cat(1,cellstr(' '),subjects(:));
            xlsout = cat(2,subjects,xlsout);
            lab_write_xls(fullfile(settings.searchfolder,['Clustering_' Foldername '_' Vars{V} '.xlsx']),xlsout);
        end
    end
end

save(fullfile(settings.searchfolder,'KmeansClustering.mat'),'Output');