function [DATA,settings,Mode] = lab_plot_read_datanew(settings,labels,Mode,defaultsel)

if isempty(Mode) | strcmp(Mode,'Select')
    disp('Labels, Matrix or File')
    Mode = questdlg('Labels, Matrix or File?','Select Input','Labels','Matrix','File','File');
elseif strcmp(Mode,'Add')
    disp('Matrix or File')
    Mode = questdlg('Matrix or File?','Select Input','Matrix','File','File');
end
if isempty(Mode)
    DATA = [];
    return
end
DATA = struct('data',[],'name','','measure','','subject','','connections',false','nodes',false,'labels',false);
if strcmp(Mode,'File')
    [DATA_file,DATA_filepath]=uigetfile('*.*','Select BrainWave or xls-file');
    if ~isempty(DATA_file) & DATA_file ~= 0
        cd(DATA_filepath);
        [~,~,Format] = lab_filename(DATA_file);
        if strcmp(Format,'xls') | strcmp(Format,'xlsx')
            try
                % try to read mappings-file
                mappings = lab_read_mappings(fullfile(DATA_filepath,DATA_file));
                if ~isfield(mappings,'mappings')
                    error('wrong file')
                end
                if mappings.mappingsChannels ~= size(labels,1)
                    disp('Abort: input files not matching')
                    DATA = [];
                    return
                end
                mapvar = questdlg('Plot FullNames, ShortNames, Single?','Selection','FullNames','ShortNames','Single','Single');
                if strcmp(mapvar,'Single')
                    for i = 1:size(mappings.mappings,2)
                        DATA(1,i).data = zeros(1,mappings.mappingsChannels);
                        DATA(1,i).data(mappings.mappings{1,i}) = 1;
                        if strcmp(mapvar,'FullNames')
                            DATA(1,i).name = mappings.mappingstitle{i};
                        else
                            DATA(1,i).name = mappings.mappingstitleS{i};
                        end
                        DATA(1,i).labels = true;
                        DATA(1,i).subject = DATA(1,i).name;
                        DATA(1,i).measure = 'Mapping';
                    end
                    settings.mappings = [];
                else
                    for i = 1:size(mappings.mappings,2)
                        DATA(1,i).data = zeros(1,size(mappings.mappings,2));
                        DATA(1,i).data(i) = 1;
                        DATA(1,i).labels = true;
                        if strcmp(mapvar,'FullNames')
                            DATA(1,i).name = mappings.mappingstitle{i};
                        else
                            DATA(1,i).name = mappings.mappingstitleS{i};
                        end
                        DATA(1,i).subject = DATA(1,i).name;
                        DATA(1,i).measure = 'Mapping';
                    end
                    settings.mappings = mappings;
                end
                settings.clustervars = mappings.mappingsChannels;
                clearvars mappings
                Mode = 'Mapping';
            catch %#ok<CTCH>
                inputdata = lab_read_xls(fullfile(DATA_filepath,DATA_file));
                inputdata = lab_correctheader(inputdata);
                try
                    tmp = find(cellfun(@isempty,inputdata));
                    if ~isempty(tmp)
                        for i = 1:length(tmp)
                            inputdata{tmp(i)} = NaN;
                        end
                    end
                    data = cell2mat(inputdata(2:end,2:end));
                catch %#ok<CTCH>
                    DATA = [];
                    return
                end
                disp('Subjects in row or column')
                Mtranspose = questdlg('Are data variables in columns or rows?',' data in columns/rows','Cancel','Rows','Columns','Columns');
                if isempty(Mtranspose) | strcmp(Mtranspose,'Cancel')
                    DATA = [];
                    return
                end
                if strcmp(Mtranspose,'Rows')
                    data = data';
                    inputdata = inputdata';
                end
                % find number of channels/locations
                if ~isfield(settings,'clustervars') | isempty(settings.clustervars)
                    [~,settings] = lab_getstructure(inputdata,settings);
                end
                if settings.clustervars == 1
                    settings.clustervars = size(data,1);
                end
                
                subjects = inputdata(1,2:end);
                NumM = floor(size(data,1)/settings.clustervars);
                for i = 1:NumM
                    tmp = strfind(inputdata{(i-1)*settings.clustervars+2,1},'_');
                    if ~isempty(tmp)
                        vars{i} = inputdata{(i-1)*settings.clustervars+2,1}(1:tmp(end)-1); %#ok<AGROW>
                    elseif settings.clustervars > 1 & NumM == 1
                        vars{i} = 'Measure'; %#ok<AGROW>
                    else
                        vars{i} = inputdata{(i-1)*settings.clustervars+2,1}; %#ok<AGROW>
                    end
                    clearvars tmp
                end
                for i = 1:size(data,2)
                    for j = 1:NumM
                        DATA(j,i).data = data((j-1)*settings.clustervars+1:j*settings.clustervars,i)';
                        if size(data,2) > 1 & ~strcmp(vars{j},'Measure')
                            DATA(j,i).name = [subjects{i} '_' vars{j}];
                        elseif size(data,2) > 1
                            DATA(j,i).name = subjects{i};
                        else
                            DATA(j,i).name = vars{j};
                        end
                        DATA(j,i).subject = subjects{i};
                        DATA(j,i).measure = vars{j};
                    end
                end
            end
        else
            [~,Result] = lab_import_bw_results(fullfile(DATA_filepath,DATA_file),[],'Single','all');
            if ~isfield(Result,'data') | isempty(Result.data)
                DATA = [];
                return
            end
            subjects = Result.subjects;
            vars = Result.results;
            for i = 1:size(Result.data,2)
                for j = 1:size(Result.data,3)
                    DATA(j,i).data = Result.data(:,i,j)';
                    if size(Result.data,2) > 1
                        DATA(j,i).name = [subjects{i} '_' vars{j}];
                    else
                        DATA(j,i).name = vars{j};
                    end
                    DATA(1,i).subject = subjects{i};
                    DATA(1,i).measure = vars{j};
                end
            end
        end
        [~,~,~,DATA_file] = lab_filename(DATA_file);
        settings.DATA_file = fullfile(DATA_filepath,DATA_file);
    else
        DATA = [];
        return
    end
elseif strcmp(Mode,'Labels')
    if isempty(labels)
        DATA = [];
        return
    end
    if ~exist('defaultsel','var')
        defaultsel = 1:length(labels);
    end
    [selection] = listdlg('PromptString','Select Regions','SelectionMode','multiple','ListString',labels,'InitialValue',defaultsel);
    if ~isempty(selection)
        for i = 1:length(selection)
            DATA(i,1).data = zeros(1,length(labels));
            DATA(i,1).data(selection(i)) = 1;
            DATA(i,1).labels = true;
            DATA(i,1).name = labels{selection(i)};
            DATA(i,1).subject = 'Labels';
            DATA(i,1).measure = labels{selection(i)};
        end
    else
        DATA = [];
        return
    end
    clearvars strlist selection
    settings.DATA_file = '';
elseif strcmp(Mode,'Matrix')
    [matrix,header,settings2] = lab_read_matrix;
    if isempty(matrix)
        DATA = [];
        return
    end
    if isfield(header,'channels')
        labels = cellstr(header.channels);
    else
        labels = cell(size(matrix,1),1);
        for i = 1:size(matrix,1)
            labels{i,1} = ['Ch' num2str(i,'%03d')];
        end
    end
    data = lab_extract_tril_wodiag(matrix,labels);
    if isfield(header,'subjects') & length(header.subjects) == size(data,2)
        subjects = header.subjects;
    else
        subjects = cell(1,size(data,2));
        for i = 1:size(data,2)
            subjects{1,i} = ['Trial' num2str(i,'%02d')];
        end
    end
    DATA_file = settings2.EEG_file;
    DATA_filepath = settings2.EEG_filepath;
    clearvars settings2
    for i = 1:size(data,2)
        DATA(1,i).data = data(:,i)';
        DATA(1,i).name = subjects{i};
        DATA(1,i).subject = subjects{i};
        DATA(1,i).measure = 'Connections';
        DATA(1,i).connections = true;
        tmp = diag(matrix(:,:,i));
        if length(unique(tmp)) > 1 | (unique(tmp) ~= 1 & unique(tmp) ~= 0)
            DATA(2,i).data = tmp(:)';
            DATA(2,i).name = subjects{i};
            DATA(2,i).subject = subjects{i};
            DATA(2,i).measure = 'Nodes';
            DATA(2,1).nodes = true;
        end
    end
    [~,~,~,DATA_file] = lab_filename(DATA_file);
    settings.DATA_file = fullfile(DATA_filepath,DATA_file);
else
    DATA = [];
    return
end