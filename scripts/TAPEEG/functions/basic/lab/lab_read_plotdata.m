% Read xls or BrainWave data for plotting
%
% [data,matrixout,cfg,skipprocessing] = lab_read_plotdata(data,labels,cfg,selectfile)
%
% written by F. Hatz 2013

function [data,matrixout,cfg,skipprocessing] = lab_read_plotdata(data,labels,cfg,selectfile)

skipprocessing = 0;

if ~exist('data','var')
    data = [];
end
matrixout = [];
header = [];
if ~exist('labels','var') | isempty(labels)
    if ~isempty(data)
        defanswer = {num2str(size(data,2))};
    else
        defanswer = {''};
    end
    answer = inputdlg({'Number of channels/regions'},'Number of variables',[1 50],defanswer);
    answer = str2num(answer{1}); %#ok<ST2NM>
    for i = 1:answer
        labels{i,1} = ['N_' num2str(i,'%03g')];
    end
    clearvars answer defanswer i
end

cfg.dofile = 0;
cfg.doregion = 0;
cfg.domatrix = 0;
cfg.domappings = 0;
if ~isfield(cfg,'dovol')
    cfg.dovol = 0;
end
if ~isfield(cfg,'dostore')
    cfg.dostore = 0;
end

if ~exist('data','var') | isempty(data)
    disp('Matrix, Selection or File')
    if exist('selectfile','var') & selectfile == 1
        button = 'File';
    else
        button = questdlg('Matrix, Selection or File?','Selection','Matrix','Selection','File','File');
    end
    if strcmp(button,'File')
        [cfg.data_file,cfg.data_filepath]=uigetfile('*.*','Select BrainWave or xls-file');
        [~,~,Format] = lab_filename(cfg.data_file);
        if ~isempty(cfg.data_file) & cfg.data_file ~= 0
            cd(cfg.data_filepath);
            if strcmp(Format,'xls') | strcmp(Format,'xlsx')
                try
                    % try to read mappings-file
                    mappings = lab_read_mappings(fullfile(cfg.data_filepath,cfg.data_file));
                    if ~isfield(mappings,'mappings')
                        error('wrong file')
                    end
                    if mappings.mappingsChannels ~= size(labels,1)
                        disp('Abort: input files not matching')
                        data = [];
                        skipprocessing = 1;
                        return
                    end
                    mapvar = questdlg('Plot FullNames, ShortNames, Single?','Selection','FullNames','ShortNames','Single','Single');
                    if strcmp(mapvar,'Single')
                        data = zeros(mappings.mappingsChannels,size(mappings.mappings,2));
                        for i = 1:size(mappings.mappings,2)
                            data(mappings.mappings{1,i},i) = 1;
                        end
                        varnames2 = mappings.mappingstitle;
                    else
                        data = eye(size(mappings.mappings,2));
                        if strcmp(mapvar,'FullNames')
                            cfg.shortnames = false;
                        else
                            cfg.shortnames = true;
                        end
                        cfg.mappings = mappings;
                        cfg.varnames = mappings.mappingstitle;
                        cfg.domappings = 1;
                    end
                    header = [];
                    cfg.doregion = 1;
                    clustervars = mappings.mappingsChannels;
                    clearvars mappings
                catch %#ok<CTCH>
                    header = lab_read_xls(fullfile(cfg.data_filepath,cfg.data_file));
                    header = lab_correctheader(header);
                    try 
                        tmp = find(cellfun(@isempty,header));
                        if ~isempty(tmp)
                            for i = 1:length(tmp)
                                header{tmp(i)} = NaN;
                            end
                        end
                        data = cell2mat(header(2:end,2:end));
                    catch %#ok<CTCH>
                        skipprocessing = 1;
                        return
                    end
                    disp('Subjects in row or column')
                    Mtranspose = questdlg('Are data variables in columns or rows?',' data in columns/rows','Cancel','Rows','Columns','Columns');
                    if isempty(Mtranspose) | strcmp(Mtranspose,'Cancel')
                        data = [];
                        skipprocessing = 1;
                        return
                    end
                    if strcmp(Mtranspose,'Rows')
                        data = data';
                        header = header';
                    end
                end
            elseif strcmp(Format,'txt') | isempty(Format)
                header = lab_import_bw_results(fullfile(cfg.data_filepath,cfg.data_file),[],'Single','all');
                try
                    data = cell2mat(header(2:end,2:end));
                catch %#ok<CTCH>
                    disp('    wrong input data')
                    data = [];
                    skipprocessing = 1;
                    return
                end
            else
                [data,header] = lab_read_data(fullfile(cfg.data_filepath,cfg.data_file));
                if isempty(data) | isempty(header) | ~isfield(header,'channels')
                    disp('    wrong input data')
                    data = [];
                    skipprocessing = 1;
                    return
                end
                data = data';
                tmp = cat(2,cellstr(header.channels),num2cell(data));
                clearvars header data
                header = {['C' num2str(size(data,1)) ' R0']};
                for j = 1:size(data,2)
                    header{1,end+1} = ['Trial' num2str(j)]; %#ok<AGROW>
                end
                header = cat(1,header,tmp);
                clearvars tmp
            end
            % find number of channels/locations
            if ~isempty(header) & size(header,2) == size(data,2)+1 & size(header,1) == size(data,1)+1
                tmp = header{1,1};
                if ~isempty(tmp) & ischar(tmp)
                    tmp = textscan(tmp,'%s');
                    tmp = tmp{1,1};
                    if strcmp(tmp{1,1}(1),'C')
                        clustervars = str2num(tmp{1,1}(2:end)); %#ok<ST2NM>
                    end
                    if size(tmp,1) > 1 & strcmp(tmp{2,1}(1),'R')
                        numresults = str2num(tmp{2,1}(2:end)); %#ok<NASGU,ST2NM>
                    end
                end
                clearvars tmp
                if ~exist('clustervars','var') | isempty(clustervars)
                    for i = 2:size(header,1)
                        tmp = strfind(header{i,1},'_');
                        if ~isempty(tmp)
                            tmp2{i-1,1} = header{i,1}(1:tmp(end)-1); %#ok<AGROW>
                        else
                            tmp2{i-1,1} = header{i,1}; %#ok<AGROW>
                        end
                    end
                    if exist('tmp2','var')
                        [~,m,~]=unique(tmp2,'last');
                        m = sort(m);
                        if length(m) == 1
                            clustervars = m(1);
                        elseif length(m) == size(tmp2,1)
                            clustervars = size(tmp2,1);
                        else
                            clustervars = m(1);
                        end
                    else
                        clustervars = [];
                    end
                    clearvars m n tmp tmp2 answer i
                end
            else
                header = [];
                if ~exist('clustervars','var')
                    clustervars = size(data,1);
                end
            end
            if clustervars ~= size(labels,1)
                answer=inputdlg({'Number of channels/locations'},'Number of channels',[1 50],{num2str(clustervars)});
                if isempty(answer)
                    skipprocessing = 1;
                    return
                end
                clustervars = str2num(answer{1,1}); %#ok<ST2NM>
            end
            if clustervars ~= size(labels,1)
                for i = 1:clustervars
                    tmp = strfind(header{i+1,1},'_');
                    if ~isempty(tmp)
                        tmp2{i,1} = header{i+1,1}(tmp(end)+1:end); %#ok<AGROW>
                    else
                        tmp2{i,1} = header{i+1,1}; %#ok<AGROW>
                    end
                end
                labels = tmp2(1:clustervars,1);
                clearvars tmp tmp2
            end
            cfg.dofile = 1;
        else
            skipprocessing = 1;
            return
        end
    elseif strcmp(button,'Selection')
        strlist = labels;
        [selection] = listdlg('PromptString','Select Regions','SelectionMode','multiple','ListString',strlist);
        if ~isempty(selection)
            data = zeros(size(selection,2),size(labels,1));
            for datanr = 1:size(selection,2)
                data(datanr,selection(datanr)) = 1;
            end
            varnames2 = labels(selection);
            cfg.doregion = 1;
        else
            data = zeros(1,size(labels,1));
        end
        clearvars strlist selection
    elseif strcmp(button,'Matrix')
        [matrix,header,cfg] = lab_read_matrix([],cfg);
        cfg.data_file = cfg.EEG_file;
        cfg.data_filepath = cfg.EEG_filepath;
        if size(matrix,3) > 1
            matrix = mean(matrix,3);
            cfg.data_file = 'Average_matrix.txt';
        end
        if isempty(matrix) | ~isnumeric(matrix) | ~size(matrix,2) == size(matrix,1)
            disp('    Abort, no valid matrix information')
            skipprocessing = 1;
            return
        else
            cfg.domatrix = 1;
            data = zeros(1,size(matrix,1));
            labels = cellstr(num2str((1:size(matrix,1))'));
        end
    else
        skipprocessing = 1;
        return
    end
elseif iscell(data)
    header = data;
    try 
        data = cell2mat(header(2:end,2:end));
        cfg.dofile = 1;
    catch %#ok<CTCH>
        disp('    wrong input data')
        data = [];
        skipprocessing = 1;
        return
    end
elseif size(data,1) == size(data,2)
    matrix = data;
    data = zeros(1,size(data,1));
    cfg.domatrix = 1;
    labels = cellstr(num2str((1:size(matrix,1))'));
elseif size(data,2) == size(labels,1)
    cfg.dofile = 1;
elseif size(data,1) == size(labels,1)
    data = data';
    cfg.dofile = 1;
end

% Control if data has to be transposed
if size(data,1) >= size(labels,1) & mod(size(data,1),size(labels,1)) < 10
    header = header';
    data = data';
end
if size(data,2) < size(labels,1) & cfg.domappings ~= 1
    disp('    Abort, mismatch data')
    skipprocessing = 1;
    return
end

if cfg.domatrix~=1 & cfg.domappings ~= 1;
    % select measures if results of different measures follow in one line
    measures = floor(size(data,2) / size(labels,1));
    for i = 1:measures
        if exist('header','var') & ~isempty(header)
            if size(header,2) > size(data,2)
                tmp = strfind(header{1,(i-1)*size(labels,1)+2},'_');
                if ~isempty(tmp)
                    strlist{i} = header{1,(i-1)*size(labels,1)+2}(1:tmp(end)-1);
                else
                    strlist{i} = header{1,(i-1)*size(labels,1)+2};
                end
            else
                tmp = strfind(header{1,(i-1)*size(labels,1)+1},'_');
                if ~isempty(tmp)
                    strlist{i} = header{1,(i-1)*size(labels,1)+1}(1:tmp(end)-1);
                else
                    strlist{i} = header{1,(i-1)*size(labels,1)+1};
                end
            end
            clearvars tmp
        else
            strlist{i} = ['Measure_' num2str(i)];
        end
    end
    if length(strlist) > 1
        [selection] = listdlg('PromptString','Select Measures','SelectionMode','multiple','ListString',strlist);
    else
        selection = 1;
    end
    if ~isempty(selection)
        datatmp = [];
        for i = 1:size(selection,2)
            datatmp(:,:,i) = data(:,(selection(i)-1)*size(labels,1)+1:selection(i)*size(labels,1)); %#ok<AGROW>
        end
        data = permute(datatmp,[3 2 1]);
        cfg.varnames = strlist(selection);
        clearvars datatmp
    else
        skipprocessing = 1;
        return
    end
    clearvars measures strlist i selection
end

if cfg.dofile == 1 & cfg.doregion ~= 1 & cfg.domatrix ~= 1 & size(data,3) > 1
    % Select measures if multiple subjects/rows
    if exist('header','var') & ~isempty(header)
        if strcmp(header{1,1},'Results') | ischar(header{2,1})
            strlist = [cellstr('--Average--') header(2:end,1)'];
        else
            strlist = [cellstr('--Average--') header(:,1)'];
        end
    else
        for i = 1:size(data,3)
            strlist{1,i} = ['Measure_' num2str(i)];
        end
        strlist = [cellstr('--Average--') strlist];
    end
    if length(strlist) > 2
        [selection] = listdlg('PromptString','Select measures/subjects','SelectionMode','multiple','ListString',strlist);
    else
        selection = 2;
    end
    if isempty(selection)
        skipprocessing = 1;
        return
    else
        selection = selection - 1;
        if min(selection) == 0 & max(selection) == 0
            data = mean(data,3);
            varnames2 = strlist(1,1);
        else
            selection = setdiff(selection,0);
            doavg = 0;
            if length(selection) > 1
                avgdialog = 0;
                for i = selection+1
                    if ~strcmp(strlist{i},'p') & ~strcmp(strlist{i},'mxp') & ~strcmp(strlist{i},'mxpV') & ...
                            ~strcmp(strlist{i},'mxpV2') & ~strcmp(strlist{i},'mxpC')
                        avgdialog = 1;
                    end
                end
                if avgdialog == 1
                    buttonavg = questdlg('Average selected measures or plot every measure?','Average measures','Cancel','Average','Single','Single');
                    if strcmp(buttonavg,'Average')
                        doavg = 1;
                    elseif ~strcmp(buttonavg,'Single')
                        skipprocessing = 1;
                        return
                    end
                end
            end
            if doavg == 1
                data = mean(data(:,:,selection),1);
                varnames2 = strlist(1,1);
            elseif doavg == 0
                data = data(:,:,selection);
                varnames2 = strlist(1,selection+1);
            end
        end
    end
    clearvars strlist selection header
end

if size(data,1) == 1 & size(data,3) > 1
    data = permute(data,[3 2 1]);
    if strcmp(cfg.varnames{1,1},'Measure_1')
        cfg.varnames = varnames2;
    else
        tmp = cfg.varnames{1,1};
        for i=1:size(varnames2,2)
            cfg.varnames{1,i} = [tmp '_' varnames2{1,i}];
        end
        clearvars tmp
    end
elseif size(data,1) == 1 & size(data,3) == 1 & exist('varnames2','var')
    if strcmp(cfg.varnames{1,1},'Measure_1')
        cfg.varnames = varnames2;
    elseif ~strcmp(varnames2{1,1},'Measure_1')
        cfg.varnames{1,1} = [varnames2{1,1} '_' cfg.varnames{1,1}];
    end
elseif size(data,1) > 1 & size(data,3) == 1 & exist('varnames2','var')
    if ~strcmp(varnames2{1,1},'Measure_1')
        for i = 1:size(data,1)
            cfg.varnames{i} = [cfg.varnames{i} '_' varnames2{1,1}];
        end
    end
elseif isfield(cfg,'varnames') & cfg.domappings ~= 1
    data = reshape(permute(data,[3 1 2]),size(data,1)*size(data,3),size(data,2));
    if exist('varnames2','var')
        varnames = [];
        for i = 1:size(cfg.varnames,2)
            for j = 1:size(varnames2,2)
                varnames = [varnames {[cfg.varnames{i} '_' varnames2{j}]}]; %#ok<AGROW>
            end
        end
        cfg.varnames = varnames;
        clearvars i j varnames
    end
end

if exist('selectfile','var') & selectfile == 1
    return
end

maxdata = size(data,1);
inverseflag = cell(1,maxdata);
for datanr = 1:maxdata
    % select settings for every datarow
    if cfg.doregion~=1 & cfg.domatrix~=1
        if strcmp(cfg.varnames{datanr},'p') | strcmp(cfg.varnames{datanr},'mxp') | ...
                strcmp(cfg.varnames{datanr},'mxpV') | strcmp(cfg.varnames{datanr},'mxpV2') | ...
                strcmp(cfg.varnames{datanr},'mxpC')
            inverseflag{datanr} = cfg.varnames{datanr};
        elseif strfind(cfg.varnames{datanr},'_p')
            inverseflag{datanr} = 'p';
        elseif strfind(cfg.varnames{datanr},'_mxp')
            inverseflag{datanr} = 'mxp';
        elseif strfind(cfg.varnames{datanr},'_mxpV')
            inverseflag{datanr} = 'mxpV';
        elseif strfind(cfg.varnames{datanr},'_mxpV2')
            inverseflag{datanr} = 'mxpV2';
        elseif strfind(cfg.varnames{datanr},'_mxpC')
            inverseflag{datanr} = 'mxpC';
        else
            inverseflag{datanr} = '';
        end
        if isfield(cfg,'varnames')
            PLOT(datanr).Name = regexprep(cfg.varnames{datanr},'_',' '); %#ok<AGROW>
        else
            PLOT(datanr).Name = ['Data ' num2str(datanr)]; %#ok<AGROW>
        end
        FORMAT{strcmp(fieldnames(PLOT),'Name')} = 'text'; %#ok<AGROW>
        if strcmp(inverseflag{datanr},'p') | strcmp(inverseflag{datanr},'mxp') | ...
                strcmp(inverseflag{datanr},'mxpV') | strcmp(inverseflag{datanr},'mxpV2') | ...
                strcmp(inverseflag{datanr},'mxpC')
            PLOT(datanr).Inverse = true; %#ok<AGROW>
        else
            PLOT(datanr).Inverse = false; %#ok<AGROW>
        end
        if PLOT(datanr).Inverse == true
            PLOT(datanr).MinValue = 0; %#ok<AGROW>
            PLOT(datanr).MaxValue = 0; %#ok<AGROW>
            PLOT(datanr).Threshold1 = 0.01; %#ok<AGROW>
            PLOT(datanr).Threshold2 = 0.05; %#ok<AGROW>
            if datanr > 1 & (strcmp(inverseflag{datanr},'mxp') | strcmp(inverseflag{datanr},'mxpV') | ...
                    strcmp(inverseflag{datanr},'mxpV2') | strcmp(inverseflag{datanr},'mxpC')) & ...
                    ~strcmp(inverseflag{datanr},inverseflag{datanr-1})
                PLOT(datanr).AddPlot = true; %#ok<AGROW>
            else
                PLOT(datanr).AddPlot = false; %#ok<AGROW>
            end
            if isfield(cfg,'LOCS')
                if PLOT(datanr).AddPlot == false
                    PLOT(datanr).Color1 = [0 1 0]; %#ok<AGROW>
                    PLOT(datanr).Color2 = [0 1 1]; %#ok<AGROW>
                else
                    PLOT(datanr).Color1 = [1 0 0]; %#ok<AGROW>
                    PLOT(datanr).Color2 = [1 0 1]; %#ok<AGROW>
                end
            else
                if PLOT(datanr).AddPlot == false
                    PLOT(datanr).Color1 = [0 0 1]; %#ok<AGROW>
                else
                    PLOT(datanr).Color1 = [1 0 0]; %#ok<AGROW>
                end
                if cfg.dovol == 1
                    PLOT(datanr).VolumePlot = false; %#ok<AGROW>
                end
            end
        elseif cfg.dofile == 1
            if (strcmp(cfg.varnames{datanr}(end),'T') | strcmp(cfg.varnames{datanr}(end),'F') |  ...
                    strcmp(cfg.varnames{datanr}(end),'Z'))
                PLOT(datanr).MinValue = -5; %#ok<AGROW>
                PLOT(datanr).MaxValue = 5; %#ok<AGROW>
                PLOT(datanr).Threshold1 = 0; %#ok<AGROW>
                PLOT(datanr).Threshold2 = 0; %#ok<AGROW>
                PLOT(datanr).Color1 = 'bluered'; %#ok<AGROW>
                PLOT(datanr).Color2 = [1 1 1]; %#ok<AGROW>
                FORMAT{strcmp(fieldnames(PLOT),'Color1')} = {'bluered','autumn','bone','colorcube','cool','copper','gray', ...
                    'hot','hsv','jet','iautumn','ibone','icolorcube','icool','icopper','igray','ihot','ihsv','ijet'}; %#ok<AGROW>
            elseif strcmp(cfg.varnames{datanr}(end),'R')
                PLOT(datanr).MinValue = -1; %#ok<AGROW>
                PLOT(datanr).MaxValue = 1; %#ok<AGROW>
                PLOT(datanr).Threshold1 = 0; %#ok<AGROW>
                PLOT(datanr).Threshold2 = 0; %#ok<AGROW>
                PLOT(datanr).Color1 = 'bluered'; %#ok<AGROW>
                PLOT(datanr).Color2 = [1 1 1]; %#ok<AGROW>
                FORMAT{strcmp(fieldnames(PLOT),'Color1')} = {'bluered','autumn','bone','colorcube','cool','copper','gray', ...
                    'hot','hsv','jet','iautumn','ibone','icolorcube','icool','icopper','igray','ihot','ihsv','ijet'}; %#ok<AGROW>
            else
                if min(data(datanr,:)) < 0
                    if min(data(datanr,:)) > -1
                        PLOT(datanr).MinValue = -max(max(abs(data(datanr,:,:)))); %#ok<AGROW>
                    else
                        PLOT(datanr).MinValue = -ceil(max(max(abs(data(datanr,:,:))))); %#ok<AGROW>
                    end
                else
                    PLOT(datanr).MinValue = min(min(data(datanr,:,:))); %#ok<AGROW>
                end
                if max(abs(data(datanr,:))) < 1
                    PLOT(datanr).MaxValue = max(max(abs(data(datanr,:,:)))); %#ok<AGROW>
                else
                    PLOT(datanr).MaxValue = ceil(max(max(abs(data(datanr,:,:))))); %#ok<AGROW>
                end
                PLOT(datanr).Threshold1 = 0; %#ok<AGROW>
                PLOT(datanr).Threshold2 = 0; %#ok<AGROW>
                if min(data(datanr,:)) < 0
                    PLOT(datanr).Color1 = 'bluered'; %#ok<AGROW>
                    PLOT(datanr).Color2 = [1 1 1]; %#ok<AGROW>
                    FORMAT{strcmp(fieldnames(PLOT),'Color1')} = {'bluered','autumn','bone','colorcube','cool','copper', ...
                        'gray','hot','hsv','jet','iautumn','ibone','icolorcube','icool','icopper','igray','ihot','ihsv','ijet'}; %#ok<AGROW>
                else
                    PLOT(datanr).ColorMode = 'color'; %#ok<AGROW>
                    FORMAT{strcmp(fieldnames(PLOT),'ColorMode')} = {'color','bluered','autumn','bone','colorcube','cool', ...
                        'copper','gray','hot','hsv','jet','iautumn','ibone','icolorcube','icool','icopper','igray','ihot','ihsv','ijet'}; %#ok<AGROW>
                    PLOT(datanr).Color1 = [1 0 0]; %#ok<AGROW>
                    PLOT(datanr).Color2 = [1 1 1]; %#ok<AGROW>
                end
            end
            PLOT(datanr).AddPlot = false; %#ok<AGROW>
            if cfg.dovol == 1
                PLOT(datanr).VolumePlot = false; %#ok<AGROW>
            end
        end
    elseif cfg.doregion == 1
        if isfield(cfg,'varnames')
            PLOT(datanr).Name = regexprep(cfg.varnames{datanr},'_',' '); %#ok<AGROW>
        else
            PLOT(datanr).Name = ['Data ' num2str(datanr)]; %#ok<AGROW>
        end
        FORMAT{strcmp(fieldnames(PLOT),'Name')} = 'text'; %#ok<AGROW>
        PLOT(datanr).Channels = find(data(datanr,:)==1); %#ok<AGROW>
        FORMAT{strcmp(fieldnames(PLOT),'Channels')} = 'vector'; %#ok<AGROW>
        if datanr == 1
            PLOT(datanr).Color1 = [0 1 0]; %#ok<AGROW>
        elseif datanr == 2
            PLOT(datanr).Color1 = [0 1 1]; %#ok<AGROW>
        elseif datanr == 3
            PLOT(datanr).Color1 = [1 0 1]; %#ok<AGROW>
        elseif datanr == 4
            PLOT(datanr).Color1 = [1 1 0]; %#ok<AGROW>
        else
            PLOT(datanr).Color1 = rand(1,3); %#ok<AGROW>
        end
        if isfield(cfg,'LOCS') & isfield(cfg.LOCS,'x')
            PLOT(datanr).MaxDistance = 0; %#ok<AGROW>
            PLOT(datanr).WeightMaxDistance = 100; %#ok<AGROW>
        end
        if datanr > 1
            PLOT(datanr).AddPlot = true; %#ok<AGROW>
        else
            PLOT(datanr).AddPlot = false; %#ok<AGROW>
        end
        if cfg.dovol == 1
            PLOT(datanr).VolumePlot = false; %#ok<AGROW>
        end
    end
end

if exist('PLOT','var')
    numdiags = ceil(size(PLOT,2) / 20);
    if ~exist('FORMAT','var')
        FORMAT = {};
    end
    for i = 1:numdiags
        Pstart = (i-1)*20+1;
        Pstop = i*20;
        if Pstop > size(PLOT,2)
            Pstop = size(PLOT,2);
        end
        [answer,skipprocessing] = inputsdlg(PLOT(Pstart:Pstop),'Plot settings',FORMAT);
        if skipprocessing == 1
            return
        else
            PLOT(Pstart:Pstop) = answer;
            clearvars answer
        end
    end
    for datanr = 1:maxdata
        if isfield(PLOT(datanr),'Inverse') & PLOT(datanr).Inverse == true
            for i = 1:size(data,3)
                if ~isempty(PLOT(datanr).Threshold1)
                    tmp = (data(datanr,:,i) - PLOT(datanr).Threshold1);
                    tmp(tmp<0) = 0;
                    if ~isempty(PLOT(datanr).Threshold2)
                        tmp = tmp / (PLOT(datanr).Threshold2 - PLOT(datanr).Threshold1);
                    end
                    tmp = tmp/2;
                    tmp(tmp>0.5) = 1;
                    data(datanr,:,i) = 1 - tmp;
                    clearvars tmp
                    PLOT(datanr).MinValue = 0; %#ok<AGROW>
                    PLOT(datanr).MaxValue = 1; %#ok<AGROW>
                elseif max(data(datanr,:,i)) <= 1 & min(data(datanr,:,i)) >= 0
                    data(datanr,:,i) = 1 - data(datanr,:,i);
                    PLOT(datanr).MinValue = 0; %#ok<AGROW>
                    PLOT(datanr).MaxValue = 1; %#ok<AGROW>
                elseif min(data(datanr,:,i)) < 0
                    data(datanr,:,i) = - data(datanr,:,i);
                    tmp = PLOT(datanr).MaxValue;
                    PLOT(datanr).MaxValue = -PLOT(datanr).MinValue; %#ok<AGROW>
                    PLOT(datanr).MinValue = -tmp; %#ok<AGROW>
                    clearvars tmp
                else
                    data(datanr,:,i) = data(datanr,:,i).^-1;
                    PLOT(datanr).MinValue = 0; %#ok<AGROW>
                    PLOT(datanr).MaxValue = round(max(data(datanr,:,i))); %#ok<AGROW>
                end
            end
        end
        if isfield(PLOT(datanr),'Channels') & ~isempty(PLOT(datanr).Channels) & isfield(PLOT(datanr),'MaxDistance')
            channels = PLOT(datanr).Channels;
            if PLOT(datanr).MaxDistance > 0
                distances = zeros(size(cfg.LOCS.x,2),size(cfg.LOCS.x,2));
                for i = 1:size(cfg.LOCS.x,2);
                    for j = 1:size(cfg.LOCS.x,2);
                        distances(i,j) = ((cfg.LOCS.x(1,i) - cfg.LOCS.x(1,j))^2 + ...
                            (cfg.LOCS.y(1,i) - cfg.LOCS.y(1,j))^2 + ...
                            (cfg.LOCS.z(1,i) - cfg.LOCS.z(1,j))^2)^0.5;
                    end
                end
                clearvars i j
                distances = distances / min(distances(distances > 0));
                sigma = (PLOT(datanr).MaxDistance^2/(-2 * log(PLOT(datanr).WeightMaxDistance/100)))^0.5;
                row = 0;
                for i = channels
                    row = row+1;
                    for j = 1:size(cfg.LOCS.x,2);
                        if distances(i,j) <= PLOT(datanr).MaxDistance
                            weights(row,j) = exp(-distances(i,j)^2 / (sigma^2 * 2)); %#ok<AGROW>
                        else
                            weights(row,j) = 0; %#ok<AGROW>
                        end
                    end
                end
                data(datanr,:) = max(weights,[],1);
                cfg.varnames{datanr} = [cfg.varnames{datanr} '_s' num2str(round(sigma*100)/100)];
                clearvars i j distances sigma row weights
            else
                data(datanr,:) = 0;
                data(datanr,channels) = 1;
            end
            clearvars channels
        elseif isfield(PLOT(datanr),'Channels') & ~isempty(PLOT(datanr).Channels)
            data(datanr,:) = 0;
            data(datanr,PLOT(datanr).Channels) = 1;
        end
    end
end

if cfg.domatrix == 1
    data = zeros(1,size(matrix,2));
    PLOT(1).Color1 = [0 0 1];
    matrixout.matrix = matrix;
end

if cfg.dostore == 1
    PLOT(1).savepictures = false;
    Prompt = cell(0,2);
    Formats = [];
    Prompt(end+1,:) = {'Save pictures' 'savepictures'};
    Formats(end+1,1).type = 'check';
    if cfg.domatrix == 1
        if size(matrix,1) < 21 & length(unique(matrix(:))) > 2
            Prompt = cat(1,Prompt,{'Value 1 (default = all)','matrixmin1'; ...
                'Value 2 (default = 50%)','matrixmin2'; ...
                'Value 3 (default = 75%)','matrixmin3'});
            matrixtmp = lab_rm_diagonal(matrix);
            PLOT(1).matrixmin1 = 0;
            PLOT(1).matrixmin2 = median(matrixtmp(:));
            PLOT(1).matrixmin3 = median(matrixtmp(matrixtmp>PLOT(1).matrixmin2));
            clearvars matrixtmp
        elseif length(unique(matrix(:))) > 2
            Prompt = cat(1,Prompt,{'Value 1 (default = 75%)','matrixmin1'; ...
                'Value 2 (default = 87.5%)','matrixmin2'; ...
                'Value 3 (default = 96.88%)','matrixmin3'});
            matrixtmp = lab_rm_diagonal(matrix);
            PLOT(1).matrixmin1 = median(matrixtmp(:));
            PLOT(1).matrixmin1 = median(matrixtmp(matrixtmp>PLOT(1).matrixmin1));
            PLOT(1).matrixmin2 = median(matrixtmp(matrixtmp>PLOT(1).matrixmin1));
            PLOT(1).matrixmin3 = median(matrixtmp(matrixtmp>PLOT(1).matrixmin2));
            PLOT(1).matrixmin3 = median(matrixtmp(matrixtmp>PLOT(1).matrixmin3));
            clearvars matrixtmp
        else
            Prompt(end+1,:) = {'Value 1 (default = 0)','matrixmin1'};
            PLOT(1).matrixmin1 = 0;
            PLOT(1).matrixmin2 = [];
            PLOT(1).matrixmin3 = [];
        end
        if isnan(PLOT(1).matrixmin1)
            PLOT(1).matrixmin1 = 0;
        end
        if isnan(PLOT(1).matrixmin2)
            PLOT(1).matrixmin2 = [];
        end
        if isnan(PLOT(1).matrixmin3)
            PLOT(1).matrixmin3 = [];
        end
        for i = 1:size(Prompt,1)-1
            Formats(end+1,1).type = 'edit'; %#ok<AGROW>
            Formats(end,1).format = 'float';
            Formats(end,1).limits = [0 inf];
        end
    end
    [PLOT(1),Cancelled] = inputsdlg(Prompt,'Save pictures with thresholds',Formats,PLOT(1));
    if Cancelled == 1
        skipprocessing = 1;
        return
    end
    if cfg.domatrix == 1
        if ~isempty(PLOT(1).matrixmin1)
            matrixout.matrix1 = matrix;
            matrixout.matrix1(matrix < PLOT(1).matrixmin1) = 0;
            PLOT(1).matrixmin1 = num2str(round(PLOT(1).matrixmin1*1000)/1000);
        end
        if ~isempty(PLOT(1).matrixmin2)
            matrixout.matrix2 = matrix;
            matrixout.matrix2(matrix < PLOT(1).matrixmin2) = 0;
            PLOT(1).matrixmin2 = num2str(round(PLOT(1).matrixmin2*1000)/1000);
        end
        if ~isempty(PLOT(1).matrixmin3)
            matrixout.matrix3 = matrix;
            matrixout.matrix3(matrix < PLOT(1).matrixmin3) = 0;
            PLOT(1).matrixmin3 = num2str(round(PLOT(1).matrixmin3*1000)/1000);
        end
    end
end
cfg.PLOT = PLOT;

end