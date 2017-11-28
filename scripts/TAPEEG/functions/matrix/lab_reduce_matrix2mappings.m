function [Result,cfg] = lab_reduce_matrix2mappings(matrixAll,cfg,nowriting)

if ~exist('nowriting','var')
    nowriting = false;
end

Result = [];
if ~exist('cfg','var')
    cfg = [];
end
if ~exist('matrixAll','var')
    [cfg.matrix_file,cfg.matrix_filepath] = uigetfile('*.xls;*.xlsx;*.txt;*.*','Select matrix file');
    if ~cfg.matrix_file == 0 | isempty(cfg.matrix_file)
        if strcmp(cfg.matrix_file(end-2:end),'txt')
            matrixAll = lab_read_txt(fullfile(cfg.matrix_filepath,cfg.matrix_file));
        elseif strcmp(cfg.matrix_file(end-2:end),'xls') | strcmp(cfg.matrix_file(end-3:end),'xlsx')
            if ispc
                matrixAll = xlsread(fullfile(cfg.matrix_filepath,cfg.matrix_file));
            else
                matrixAll = xlsread(fullfile(cfg.matrix_filepath,cfg.matrix_file),1,'','basic');
            end
        else
            try
                matrix = lab_read_bw_matrices(fullfile(cfg.matrix_filepath,cfg.matrix_file));
                filenames = matrix.name;
                matrixAll = matrix.matrix;
                clearvars matrix
            catch %#ok<CTCH>
                matrixAll = [];
            end
        end
    else
        return
    end
elseif isfield(cfg,'Output_filepath')
    cfg.matrix_file = cfg.Output_file;
    cfg.matrix_filepath = cfg.Output_filepath;
end

if isempty(matrixAll) | size(matrixAll,1) ~= size(matrixAll,2)
    disp('    Abort, no valid matrix file')
    return
end

if ~isfield(cfg,'MATRIX') | ~isfield(cfg.MATRIX,'Mappings') | isempty(cfg.MATRIX.Mappings)
    cfg.EXTRA.numdatachans = size(matrixAll,2);
    cfg.MATRIX.Mappings = lab_load_mappings([],cfg);
end
mappings = cfg.MATRIX.Mappings;
if ~isfield(cfg.MATRIX,'MappingsMode') | isempty(cfg.MATRIX.MappingsMode)
    cfg.MATRIX.MappingsMode = 'Average';
end

if mappings.mappingsChannels ~= size(matrixAll,1)
    disp('   Skip Mappings, number of channels not matching')
    return
end

if isfield(cfg,'matrix_filepath') & size(matrixAll,3) > 1
    warning off %#ok<WNOFF>
    mkdir(fullfile(cfg.matrix_filepath,'Matrix2mappings'));
    warning off %#ok<WNOFF>
    cfg.Output_filepath = fullfile(cfg.matrix_filepath,'Matrix2mappings');
elseif isfield(cfg,'matrix_filepath')
    cfg.Output_filepath = cfg.matrix_filepath;
end

disp(['   Reduce Matrix to Mappings (number: ' num2str(size(mappings.mappings,2)) ')'])

for nmatrix = 1:size(matrixAll,3)
    matrix = matrixAll(:,:,nmatrix);
    if exist('filenames','var')
        cfg.matrix_file = filenames{1,nmatrix};
    end
    if ~isempty(mappings) & isstruct(mappings)
        numout = size(mappings.mappings,2);
        if strcmp(cfg.MATRIX.MappingsMode,'Average')
            matrix(1:size(matrix,1)+1:end) = NaN;
            matrixout = zeros(numout,numout);
            for i = 1:numout
                for j = 1:numout
                    tmp = matrix(mappings.mappings{i},mappings.mappings{j});
                    tmp = tmp(~isnan(tmp));
                    matrixout(i,j) = mean(tmp(:));
                end
            end
        elseif strcmp(cfg.MATRIX.MappingsMode,'Degree')
            Degrees = degrees_und(matrix);
            Index = zeros(1,numout);
            for i = 1:numout
                Index(1,i) = find(Degrees(mappings.mappings{i}) == max(Degrees(mappings.mappings{i})),1);
            end
            matrixout = matrix(Index,Index);
            matrixout(1:numout+1:end) = 0;
            clearvars Index Degrees
        elseif strcmp(cfg.MATRIX.MappingsMode,'Betweenness')
            matrixI = lab_invmatrix(matrix);
            [~,Betweenness] = edge_betweenness_wei(matrixI);
            Index = zeros(1,numout);
            for i = 1:numout
                Index(1,i) = find(Betweenness(mappings.mappings{i}) == max(Betweenness(mappings.mappings{i})),1);
            end
            matrixout = matrix(Index,Index);
            matrixout(1:numout+1:end) = 0;
            clearvars matrixI Index Betweenness
        end
        Result.matrix(:,:,nmatrix) = matrixout;
        header.channels = char(mappings.mappingstitleS);
        Result.channels = header.channels;
        redflag = true;
    else
        matrixout = matrix;
        redflag = false;
    end
    
    if isfield(cfg,'matrix_file')
        [~,~,~,cfg.Output_file] = lab_filename(cfg.matrix_file);
        if redflag == true
            cfg.Output_file = [cfg.Output_file '_ToMap'];
        end
    end
    if isfield(cfg,'Output_filepath') & nowriting == false
        if strcmp(cfg.Output_file(end-5:end),'matrix')
            matrixfileout = fullfile(cfg.Output_filepath,[cfg.Output_file '.txt']);
        else
            matrixfileout = fullfile(cfg.Output_filepath,[cfg.Output_file '_matrix.txt']);
        end
        if exist(matrixfileout,'file')
            delete(matrixfileout);
        end
        dlmwrite(matrixfileout,matrixout,'delimiter','\t','precision', 6);
    end
end
