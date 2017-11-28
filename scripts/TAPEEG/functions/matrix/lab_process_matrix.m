% Processing of matrices, script called by lab_tapeeg
%
% [MatrixOut,Result,cfg] = lab_process_matrix(Matrix,header,cfg,nowriting,nograph)
%
% Written by F. Hatz 2013 Neurology Basel

function [MatrixOut,Result,cfg] = lab_process_matrix(Matrix,header,cfg,nowriting,nograph)
    
MatrixOut = [];
Result = [];

global MatrixAll

if ~exist('nograph','var')
    nograph = false;
end
if ~exist('nowriting','var')
    nowriting = false;
end
if ~exist('cfg','var')
    cfg = [];
end
if ~isfield(cfg,'lastfile')
    cfg.lastfile = true;
end

if (~isfield(cfg,'Output_filepath') | isempty(cfg.Output_filepath)) & isfield(cfg,'EEG_filepath')
    cfg.Output_filepath = cfg.EEG_filepath;
    cfg.Output_file = cfg.EEG_file;
end

if isfield(cfg,'Output_file')
    if exist('header','var') & isfield(header,'subjects') & ~isempty(header.subjects)
        Filenames = header.subjects(:)';
        Patients = header.subjects(:)';
    elseif isfield(cfg.MATRIX,'subjectname')
        patient = lab_subjectname(cfg.Output_file,cfg.MATRIX);
        if size(Matrix,3) > 1
            Filenames = cell(1,size(Matrix,3));
            for i = 1:size(Matrix,3)
                Filenames{1,i} = [patient '_' num2str(i) '.txt'];
            end
        else
            Filenames{1} = [patient '.txt'];
        end
        Patients = repmat(cellstr(patient),1,size(Matrix,3));
        clearvars patient i
    elseif isfield(cfg,'patient')
        patient = cfg.patient;
        if size(Matrix,3) > 1
            Filenames = cell(1,size(Matrix,3));
            for i = 1:size(Matrix,3)
                Filenames{1,i} = [patient '_' num2str(i) '.txt'];
            end
            clearvars i
        else
            Filenames{1} = [patient '.txt'];
        end
        Patients = repmat(cellstr(cfg.patient),1,size(Matrix,3));
        clearvars patient
    else
        patient = cfg.Output_file;
        if size(Matrix,3) > 1
            Filenames = cell(1,size(Matrix,3));
            for i = 1:size(Matrix,3)
                Filenames{1,i} = [patient '_' num2str(i) '.txt'];
            end
            clearvars i
        else
            Filenames{1} = [patient '.txt'];
        end
        Patients = repmat(cellstr(cfg.Output_file),1,size(Matrix,3));
        clearvars patient
    end
    if isfield(header,'freqband') & ~isempty(header.freqband)
        cfg.Output_band = ['F' num2str(header.freqband(1)) 'F' num2str(header.freqband(2))];
    end
else
    Filenames = {};
    Patients = {};
end

if isfield(cfg,'Output_filepath') & ~isempty(cfg.Output_filepath) & nowriting == false
    Output_filepath = cfg.Output_filepath;
    Output_fileB = cfg.Output_file;
    Output_filepathB = cfg.Output_filepath;
end

if ~isfield(cfg,'MATRIX')
    if ~isempty(Filenames)
        cfg.Output_file = Filenames{1};
    end
    cfg.EXTRA.numdatachans = size(Matrix,2);
    if size(Matrix,3) > 1
        [cfg,skipprocessing] = lab_set_process_matrix(cfg,false,false,false,true,header);
    else
        [cfg,skipprocessing] = lab_set_process_matrix(cfg,false,false,false,false,header);
    end
    if skipprocessing == 1
        Result = [];
        MatrixOut = [];
        return
    end
end

if exist('Output_filepath','var') & ~isempty(Output_filepath) & nowriting == false
    if isfield(cfg.MATRIX,'outputfolder') & ~isempty(cfg.MATRIX.outputfolder)
        warning off %#ok<WNOFF>
        mkdir(fullfile(Output_filepath,cfg.MATRIX.outputfolder));
        warning on %#ok<WNON>
        Output_filepath = fullfile(Output_filepath,cfg.MATRIX.outputfolder);
    end
end

if isfield(cfg.MATRIX,'nummatrices') & ~isempty(cfg.MATRIX.nummatrices) & size(Matrix,3) > cfg.MATRIX.nummatrices
    Matrix = Matrix(:,:,1:cfg.MATRIX.nummatrices);
end

if isfield(cfg.MATRIX,'exclude') & isnumeric(cfg.MATRIX.exclude) & ~isempty(cfg.MATRIX.exclude)
    include = setdiff(1:size(Matrix,1),cfg.MATRIX.exclude);
    if ~isempty(include) & length(include) < size(Matrix,1);
        disp(['   Reduce to ' num2str(length(include)) ' channels'])
        Matrix = Matrix(include,include,:);
    end
end

if isfield(cfg.MATRIX,'donormalize') & ischar(cfg.MATRIX.donormalize) & strcmp(cfg.MATRIX.donormalize,'over all')
    disp(['   Normalize Matrix (' num2str(size(Matrix,3)) ' matrices)'])
    Matrix = lab_normalize_matrix(Matrix);
end

if isfield(cfg.MATRIX,'average') & ischar(cfg.MATRIX.average)
    if strcmp(cfg.MATRIX.average,'subject')
        Matrix = mean(Matrix,3);
    elseif strcmp(cfg.MATRIX.average,'folder') & ~isfield(cfg,'lastfilefolder')
        disp('  Average by folder not possible, replace by average all');
        cfg.MATRIX.average = 'all';
    end
end

nowriting_flag = nowriting;
for Nmatrix = 1:size(Matrix,3)
    skipprocessing = 0;
    nowriting = nowriting_flag;
    matrix = Matrix(:,:,Nmatrix);
    
    if ~isempty(Filenames)
        Filename = Filenames{Nmatrix};
        [~,~,~,cfg.Output_fileS] = lab_filename(Filename);
        if isfield(cfg,'Output_band')
            cfg.Output_file = [cfg.Output_fileS '_' cfg.Output_band '.sef'];
            cfg.Output_fileS = [cfg.Output_fileS '_' cfg.Output_band];
        else
            cfg.Output_file = [cfg.Output_fileS '.sef'];
        end
        cfg.patient = Patients{Nmatrix};
        disp(['   Process Matrix: ' cfg.Output_fileS])
    end
    
    if isfield(cfg.MATRIX,'donormalize') & ischar(cfg.MATRIX.donormalize) & strcmp(cfg.MATRIX.donormalize,'single')
        if isfield(cfg,'Output_fileS')
            disp(['   Normalize ' cfg.Output_fileS])
        end
        matrix = lab_normalize_matrix(matrix);
    end
    
    if isfield(cfg.MATRIX,'dosymmetrical') & cfg.MATRIX.dosymmetrical == true
        if isfield(cfg,'Output_fileS')
            disp(['   Transform ' cfg.Output_fileS ' to symmetric matrix'])
        end
        matrix = (matrix + tril(matrix)' + triu(matrix)') ./ 2;
        if isfield(cfg,'Output_file')
            cfg.Output_fileS = [cfg.Output_fileS '_SYM'];
            cfg.Output_file = [cfg.Output_fileS '.sef'];
        end
    end
    
    if isfield(cfg.MATRIX,'dodirectional') & cfg.MATRIX.dodirectional == true
        if isfield(cfg,'Output_fileS')
            disp(['   Transform ' cfg.Output_fileS ' to directional matrix'])
        end
        matrix = matrix ./ (matrix + matrix');
        matrix(1:size(matrix,1)+1:end) = 0;
        if isfield(cfg,'Output_file')
            cfg.Output_fileS = [cfg.Output_fileS '_DIR'];
            cfg.Output_file = [cfg.Output_fileS '.sef'];
        end
    end
    
    if isfield(cfg.MATRIX,'minzero') & cfg.MATRIX.minzero == true
        if isfield(cfg,'Output_fileS')
            disp(['   Set values < 0 of ' cfg.Output_fileS ' to zero'])
        end
        matrix(matrix < 0) = 0;
    end
    
    if isfield(cfg.MATRIX,'DIST') & ~isempty(cfg.MATRIX.DIST)
        Distance_file = '';
        DistMatrix = [];
        if isfield(cfg.MATRIX.DIST,'input') & strcmp(cfg.MATRIX.DIST.input,'off')
            if isfield(cfg.MATRIX.DIST,'matrix') & ~isempty(cfg.MATRIX.DIST.matrix)
                DistMatrix = cfg.MATRIX.DIST.matrix;
                Distance_file = 'Default';
                disp('   Correct for distances with Template')
            end
        else
            if strcmp(cfg.MATRIX.DIST.input,'file') & isfield(cfg,'EEG_file')
                [~,~,~,EEG_fileS] = lab_filename(cfg.EEG_file);
                Distance_file = fullfile(cfg.EEG_filepath,[EEG_fileS '_DistMatrix.txt']);
            elseif strcmp(cfg.MATRIX.DIST.input,'patient') & isfield(cfg,'patient')
                Distance_file = fullfile(cfg.EEG_filepath,[cfg.patient '_DistMatrix.txt']);
            end
            if ~isempty(Distance_file) & exist(Distance_file,'file');
                DistMatrix = lab_read_data(Distance_file);
                disp(['   Correct for distances with file: ' lab_filename(Distance_file)])
            end
        end
        if ~isempty(DistMatrix)
            if ~isfield(cfg.MATRIX.DIST,'threshold') | isempty(cfg.MATRIX.DIST.threshold)
                cfg.MATRIX.DIST.treshold = 5;
            end
            if ~isfield(cfg.MATRIX.DIST,'mode') | isempty(cfg.MATRIX.DIST.mode)
                cfg.MATRIX.DIST.mode = '1';
            end
            switch cfg.MATRIX.DIST.mode
                case '1'
                    disp(['   Correct distance (' Distance_file ') < ' num2str(cfg.MATRIX.DIST.threshold) ' to 1'])
                    matrix(DistMatrix <= cfg.MATRIX.DIST.threshold) = 1;
                case 'max'
                    disp(['   Correct distance (' Distance_file ') < ' num2str(cfg.MATRIX.DIST.threshold) ' to max value'])
                    matrix(DistMatrix <= cfg.MATRIX.DIST.threshold) = max(matrix(:));
            end
        else
            Distance_file = 'Disabled';
        end
    end
    
    if isfield(cfg.MATRIX,'rankmatrix') & cfg.MATRIX.rankmatrix == true
        if ~isfield(cfg.MATRIX,'rankorder') | isempty(cfg.MATRIX.rankorder)
            cfg.MATRIX.rankorder = 5;
        end
        if isfield(cfg,'Output_fileS')
            disp(['   Rank ' cfg.Output_fileS ' with order ' num2str(cfg.MATRIX.rankorder)])
        end
        matrix = lab_rank_matrix(matrix,cfg.MATRIX.rankorder);
    end
    
    if isfield(cfg.MATRIX,'average') & ischar(cfg.MATRIX.average) & strcmp(cfg.MATRIX.average,'all')
        if isfield(cfg,'Output_band') & ~isempty(cfg.Output_band)
            Output_band = cfg.Output_band;
        else
            Output_band = 'Fall';
        end
        if isfield(header,'measure') & ~isempty(header.measure)
            Measure = header.measure;
        else
            Measure = 'Measure';
        end
        if ~isfield(MatrixAll,Output_band) | ~isfield(MatrixAll.(Output_band),Measure)
            MatrixAll.(Output_band).(Measure) = [];
        end
        if isempty(MatrixAll.(Output_band).(Measure))
            if isfield(cfg,'Output_fileS')
                disp(['   Store ' cfg.Output_fileS ' for grand average matrix'])
            end
            MatrixAll.(Output_band).(Measure) = matrix;
        elseif size(MatrixAll.(Output_band).(Measure),1) == size(matrix,1) & ...
                size(MatrixAll.(Output_band).(Measure),2) == size(matrix,2)
            if isfield(cfg,'Output_fileS')
                disp(['   Store ' cfg.Output_fileS ' for grand average matrix'])
            end
            MatrixAll.(Output_band).(Measure) = cat(3,MatrixAll.(Output_band).(Measure),matrix);
        elseif isfield(cfg,'Output_fileS')
            disp(['   Skip matrix (' cfg.Output_fileS '), size not consistent with previous matrices'])
        end
        if cfg.lastfile == true & Nmatrix == size(Matrix,3)
            matrix = mean(MatrixAll.(Output_band).(Measure),3);
            if isfield(cfg.MATRIX,'dosymmetrical') & cfg.MATRIX.dosymmetrical == true
                cfg.Output_file = 'GrandAverage_SYM.txt';
            else
                cfg.Output_file = 'GrandAverage.txt';
            end
        else
            skipprocessing = 1;
        end
    elseif isfield(cfg.MATRIX,'average') & ischar(cfg.MATRIX.average) & strcmp(cfg.MATRIX.average,'folder')
        if isfield(cfg,'Output_band') & ~isempty(cfg.Output_band)
            Output_band = cfg.Output_band;
        else
            Output_band = 'Fall';
        end
        if isfield(header,'measure') & ~isempty(header.measure)
            Measure = header.measure;
        else
            Measure = 'Measure';
        end
        if ~isfield(MatrixAll,Output_band) | ~isfield(MatrixAll.(Output_band),Measure)
            MatrixAll.(Output_band).(Measure) = [];
        end
        if isempty(MatrixAll.(Output_band).(Measure))
            if isfield(cfg,'Output_fileS')
                disp(['   Store ' cfg.Output_fileS ' for folder average matrix'])
            end
            MatrixAll.(Output_band).(Measure) = matrix;
        elseif size(MatrixAll.(Output_band).(Measure),1) == size(matrix,1) & ...
                size(MatrixAll.(Output_band).(Measure),2) == size(matrix,2)
            if isfield(cfg,'Output_fileS')
                disp(['   Store ' cfg.Output_fileS ' for folder average matrix'])
            end
            MatrixAll.(Output_band).(Measure) = cat(3,MatrixAll.(Output_band).(Measure),matrix);
        elseif isfield(cfg,'Output_fileS')
            disp(['   Skip matrix (' cfg.Output_fileS '), size not consistent with previous matrices'])
        end
        if cfg.lastfilefolder == true & Nmatrix == size(Matrix,3)
            matrix = mean(MatrixAll.(Output_band).(Measure),3);
            MatrixAll.(Output_band).(Measure) = [];
            tmp.subjectname = -1;
            Name = lab_subjectname(fullfile(cfg.Output_filepath,cfg.Output_file),tmp);
            if isempty(Name)
                Name = 'FolderAverage';
            end
            if isfield(cfg.MATRIX,'dosymmetrical') & cfg.MATRIX.dosymmetrical == true
                cfg.Output_file = [Name '_SYM.txt'];
            else
                cfg.Output_file = [Name '.txt'];
            end
        else
            skipprocessing = 1;
        end
    end
    
    if skipprocessing == 0
        if isfield(cfg.MATRIX,'Mappings') & ~isempty(cfg.MATRIX.Mappings)
            if isfield(cfg,'Output_filepath') & nowriting == false
                cfg.Output_filepath = fullfile(Output_filepath,'Matrix2Map');
                warning off %#ok<WNOFF>
                mkdir(cfg.Output_filepath);
                warning on %#ok<WNON>
            end
            [Result,cfg] = lab_reduce_matrix2mappings(matrix,cfg,nowriting);
            if ~isempty(Result)
                matrix = Result.matrix;
                header.channels = Result.channels;
            end
            
            % write Mappings-Connectivities to xls
            if isfield(cfg.MATRIX,'MappingsWrite') & cfg.MATRIX.MappingsWrite == true & ...
                    isfield(cfg,'Output_filepath') & nowriting == false
                [~,~,~,FilenameS] = lab_filename(cfg.Output_file);
                if isfield(cfg,'Output_band') & ~isempty(cfg.Output_band)
                    if ~isfield(cfg,'Output_filepath') & isfield(cfg,'settings_path')
                        connectionsfileout = fullfile(cfg.settings_path,[cfg.Output_band '_Connectivities.xlsx']);
                    elseif isfield(cfg,'Output_filepath')
                        connectionsfileout = fullfile(cfg.Output_filepath,[cfg.Output_band '_Connectivities.xlsx']);
                    end
                else
                    if ~isfield(cfg,'Output_filepath') & isfield(cfg,'settings_path')
                        connectionsfileout = fullfile(cfg.settings_path,'Connectivities.xlsx');
                    elseif isfield(cfg,'Output_filepath')
                        connectionsfileout = fullfile(cfg.Output_filepath,'Connectivities.xlsx');
                    end
                end
                [connections,labels] = lab_extract_tril(matrix,cfg.MATRIX.Mappings.mappingstitleS);
                if exist(connectionsfileout,'file');
                    xlsout = lab_read_xls(connectionsfileout);
                else
                    xlsout = cat(1,{''},labels);
                end
                if size(connections,1) + 1 == size(xlsout,1)
                    xlsout{1,end+1} = FilenameS; %#ok<AGROW>
                    xlsout(2:end,end) = num2cell(connections);
                    lab_write_xls(connectionsfileout,xlsout);
                end
            end
            nowriting = true;
        end
        
        if isfield(cfg.MATRIX,'binary') & cfg.MATRIX.binary == true
            if ~isfield(cfg.MATRIX,'binarythreshold') | isempty(cfg.MATRIX.binarythreshold)
                cfg.MATRIX.binarythreshold = 0.5;
            end
            if ~isfield(cfg.MATRIX,'binarymode') | isempty(cfg.MATRIX.binarymode)
                cfg.MATRIX.binarymode = 'fixed';
            end
            matrix = lab_matrix2binary(matrix,cfg.MATRIX.binarythreshold,cfg.MATRIX.binarymode);
        end
        
        if exist('Output_filepath','var') & nowriting == false
            [~,~,~,FilenameS] = lab_filename(cfg.Output_file);
            matrixfileout = fullfile(cfg.Output_filepath,[FilenameS '_matrix.txt']);
            if exist(matrixfileout,'file')
                delete(matrixfileout);
            end
            dlmwrite(matrixfileout,matrix,'delimiter','\t','precision', 6);
            fig1 = lab_plot_matrix(matrix,true);
            lab_print_figure([matrixfileout(1:end-3) '.jpg'],fig1);
            close(fig1);
        end
        
        if ~isempty(matrix) & isfield(cfg.MATRIX,'PLOT') & ~isempty(cfg.MATRIX.PLOT)
            % Plot matrix
            if ~isempty(Filenames)
                Filename = Filenames{Nmatrix};
                [~,~,~,cfg.Output_fileS] = lab_filename(Filename);
                cfg.Output_filepath = Output_filepath;
                if isfield(cfg,'Output_band')
                    cfg.Output_file = [cfg.Output_fileS '_' cfg.Output_band '.sef'];
                    cfg.Output_fileS = [cfg.Output_fileS '_' cfg.Output_band];
                else
                    cfg.Output_file = [cfg.Output_fileS '.sef'];
                end
                [~,~,~,cfg.Output_fileS] = lab_filename(cfg.Output_file);
            end
            disp(['    Plot matrix: ' cfg.Output_file])
            settmp.handleF = figure('Visible','off','Color',[1 1 1],'Renderer','zbuffer','Menubar','none');
            settmp.PLOT_file = [];
            settmp.close = 0;
            settmp.LOCS = cfg.MATRIX.PLOT.LOCS;
            if isfield(cfg.MATRIX,'Mappings') & ~isempty(cfg.MATRIX.Mappings)
                settmp.Mappings = cfg.MATRIX.Mappings;
            else
                settmp.Mappings = cfg.MATRIX.PLOT.Mappings;
            end
            if ~isempty(settmp.Mappings)
                settmp.Mappings.shortnames = true;
            end
            PLOT = cfg.MATRIX.PLOT;
            PLOT.ColorE = lab_create_cmap(PLOT.ColorE);
            PLOT.Color = lab_create_cmap(PLOT.Color);
            if isfield(cfg.MATRIX.PLOT,'nodemethod')
                nodes = lab_plot_create_nodes(matrix,cfg.MATRIX.PLOT);
            else
                nodes = zeros(1,size(matrix,1));
            end
            matrixtmp = matrix;
            matrixtmp(1:size(matrix,1)+1:end) = nodes;
            if max(nodes(:)) == 0 & min(nodes(:)) == 0
                PLOT.MinValue = 0;
                PLOT.MaxValue = 1;
            elseif PLOT.MinValue == 0 & PLOT.MaxValue == 0
                PLOT.MinValue = min(nodes(:));
                PLOT.MaxValue = max(nodes(:));
            end
            settmp = lab_plot_chans(matrixtmp,PLOT,settmp);
            lab_print_figure(fullfile(cfg.Output_filepath,[cfg.Output_fileS '.tif']),settmp.handleF);
            close(settmp.handleF);
            clearvars PLOT settmp nodes matrixtmp
        end
        
        if exist('Output_filepath','var') & isfield(cfg,'Output_fileS') & nowriting == false
            % Write verbose file (*.vrb)
            fid=fopen(fullfile(Output_filepath,[cfg.Output_fileS '.vrb']),'w');
            fprintf(fid,'Process Matrix\n');
            fprintf(fid,['EEG-file: ' cfg.Output_fileS]);
            fprintf(fid,'\n\n');
            if isfield(cfg.MATRIX,'donormalize')
                fprintf(fid,['Normalize Matrix: ' cfg.MATRIX.donormalize '\n\n']);
            end
            if isfield(cfg.MATRIX,'dosymmetrical') & cfg.MATRIX.dosymmetrical == true
                fprintf(fid,'Do symmetrical: on\n\n');
            else
                fprintf(fid,'Do symmetrical: off\n\n');
            end
            if isfield(cfg.MATRIX,'DIST') & ~isempty(cfg.MATRIX.DIST)
                fprintf(fid,'Correct Distance:\n');
                fprintf(fid,['  File: ' lab_filename(Distance_file) '\n']);
                fprintf(fid,['  Threshold: ' num2str(cfg.MATRIX.DIST.threshold) '\n']);
                fprintf(fid,['  Value to set: ' cfg.MATRIX.DIST.mode '\n']);
                fprintf(fid,['  Matrix by name: ' cfg.MATRIX.DIST.input '\n\n']);
            end
            if isfield(cfg.MATRIX,'rankmatrix') & cfg.MATRIX.rankmatrix == true
                fprintf(fid,['Rank matrix: on (Order: ' num2str(cfg.MATRIX.rankorder) ')\n\n']);
            else
                fprintf(fid,'Rank matrix: off\n\n');
            end
            if isfield(cfg.MATRIX,'grandaverage') & cfg.MATRIX.grandaverage == true
                fprintf(fid,'Grand average: on\n\n');
            else
                fprintf(fid,'Grand average: off\n\n');
            end
            if isfield(cfg.MATRIX,'Mappings') & ~isempty(cfg.MATRIX.Mappings)
                fprintf(fid,['Mappings: on (' sprintf('%s|',cfg.MATRIX.Mappings.mappingstitle{:}) ')\n\n']);
            else
                fprintf(fid,'Mappings: off\n\n');
            end
            if isfield(cfg.MATRIX,'GRAPH') & ~isempty(cfg.MATRIX.GRAPH)
                fprintf(fid,'Graph analysis: on\n\n');
            else
                fprintf(fid,'Graph analysis: off\n\n');
            end
            fclose(fid);
        end
        
        if ~isempty(matrix) & isfield(cfg.MATRIX,'GRAPH') & ~isempty(cfg.MATRIX.GRAPH) & nograph == false
            if ~isempty(Filenames)
                timestamp = cellstr(cfg.Output_fileS);
                Filename = Filenames{Nmatrix};
                [~,~,~,cfg.Output_fileS] = lab_filename(Filename);
                cfg.Output_filepath = Output_filepath;
                if isfield(cfg,'Output_band')
                    cfg.Output_file = [cfg.Output_fileS '_' cfg.Output_band '.sef'];
                    cfg.Output_fileS = [cfg.Output_fileS '_' cfg.Output_band];
                else
                    cfg.Output_file = [cfg.Output_fileS '.sef'];
                end
                [~,~,~,cfg.Output_fileS] = lab_filename(cfg.Output_file);
            else
                timestamp = {};
            end
            if Nmatrix == size(Matrix,3)
                cfg.lastmatrix = true;
            elseif ~isempty(Patients) & ~strcmp(Patients{Nmatrix},Patients{Nmatrix+1})
                cfg.lastmatrix = true;
            else
                cfg.lastmatrix = false;
            end
            [Rtmp,cfg.MATRIX,cfg] = lab_graphanalysis(matrix,header,cfg.MATRIX,cfg,timestamp);
            cfg = rmfield(cfg,'lastmatrix');
            if ~isempty(Rtmp)
                if size(Matrix,3) > 1
                    Result.(['Graph' num2str(Nmatrix)]) = Rtmp;
                else
                    Result.Graph = Rtmp;
                end
            end
            clearvars Rtmp
        end
        
        if exist('MatrixOut','var') & ~isempty(MatrixOut) & size(MatrixOut,1) == size(matrix,1)
            MatrixOut(:,:,Nmatrix) = matrix; %#ok<AGROW>
        else
            MatrixOut = matrix;
        end
    end
end

if exist('Output_filepath','var')
    cfg.Output_filepath = Output_filepathB;
    cfg.Output_file = Output_fileB;
    [~,~,~,cfg.Output_fileS] = lab_filename(cfg.Output_file);
end

end