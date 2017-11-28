% Read xls or BrainWave data for plotting
%
% [data,cfg,skipprocessing] = lab_plot_read_data(cfg,Labels,Mode)
%
% written by F. Hatz 2014

function [DATA,cfg] = lab_plot_read_data(DATA,cfg,Labels,Mode,doIS)

global DATA_ALL
    
if ~exist('doIS','var')
    doIS = 0;
end
if ~exist('Mode','var')
    Mode = '';
end
if ~exist('Labels','var')
    Labels = [];
    numchans = [];
else
    numchans = length(Labels);
end
if ~exist('cfg','var')
    cfg = [];
end
if ~exist('DATA','var')
    DATA = [];
end

if ~isempty(DATA_ALL)
    DATA = DATA_ALL;
end

settings = [];
if ~isempty(DATA) & strcmp(Mode,'Labels')
    dataerror = false;
    for i = 1:length(DATA)
        tmp = find(strcmp(Labels,DATA(i).measure));
        if isempty(tmp)
            dataerror = true;
        else
            defaultsel(i) = tmp; %#ok<AGROW>
        end
    end
    if dataerror == true
        DATA = [];
    else
        [DATA,settings,Mode] = lab_plot_read_datanew(settings,Labels,Mode,defaultsel);
    end
elseif isempty(DATA)
    [DATA,settings,Mode] = lab_plot_read_datanew(settings,Labels,Mode);
end

if isempty(DATA)
    return
end

% Match to number of electrodes
if ~isempty(numchans)
    DATA = lab_plot_split_connections(DATA,[],numchans);
end

% Extract connections
IsConnection = find_connections(DATA);

Prompt = cell(0,2);
Formats = [];
if (size(DATA,1) > 1 & (strcmp(Mode,'File') | strcmp(Mode,'Matrix'))) | strcmp(Mode,'Add')
    Nm = strcmp(fieldnames(DATA),'measure');
    tmp = struct2cell(DATA);
    tmp = tmp(Nm,:,1);
    measures = tmp(:);
    settings.selectmeasures = 1:length(measures);
    clearvars tmp Nm
    
    Prompt(end+1,:) = {'Measures','selectmeasures'};
    Formats(end+1,1).type = 'list';
    Formats(end,1).style = 'listbox';
    Formats(end,1).items = measures;
    Formats(end,1).limits = [0 2]; % multi-select
    Formats(end,1).size = [180 250];
else
    settings.selectmeasures = 1;
    measures = cellstr(DATA(1).measure);
end
if (size(DATA,2) > 1 & (strcmp(Mode,'File') | strcmp(Mode,'Matrix'))) | strcmp(Mode,'Add')
    Nm = strcmp(fieldnames(DATA),'subject');
    tmp = struct2cell(DATA);
    tmp = tmp(Nm,1,:);
    subjects = tmp(:);
    if max(strcmp(subjects,'p')) == 1 | max(strcmp(subjects,'mxp')) == 1 | max(strcmp(subjects,'mxpV')) == 1 | ...
            max(strcmp(subjects,'mxpV2')) == 1 | max(strcmp(subjects,'mxpC')) == 1 | max(strcmp(subjects,'mxpS')) == 1
        doavg = false;
        settings.selectsubjects = 1:length(subjects);
    else
        doavg = true;
        subjects = cat(1,cellstr('--Average--'),subjects);
        settings.selectsubjects = 1;
    end
    clearvars tmp Nm
    
    Prompt(end+1,:) = {'Subjects/Trials','selectsubjects'};
    Formats(end+1,1).type = 'list';
    Formats(end,1).style = 'listbox';
    Formats(end,1).items = subjects;
    Formats(end,1).limits = [0 2]; % multi-select
    Formats(end,1).size = [180 250];
else
    settings.selectsubjects = 1;
    doavg = true;
    subjects = cellstr(DATA(1).subject); %#ok<NASGU>
end

if ~isempty(Prompt)
    Formats(end+1,1).type = 'none';
    Formats(end+1,1).type = 'none';
    
    if max(IsConnection) == true
        Prompt(end+1,:) = {'Add graph measure',''};
        Formats(end+1,1).type = 'button';
        Formats(end,1).style = 'pushbutton';
        Formats(end,1).size = [140 25];
        Formats(end,1).callback = {@add_graphmeasure,'','@selectmeasures'};
        Formats(end,1).span = [1 2];
    end
    if strcmp(Mode,'File') | strcmp(Mode,'Matrix') | strcmp(Mode,'Add')
        Prompt(end+1,:) = {'Add Data',''};
        Formats(end+1,1).type = 'button';
        Formats(end,1).style = 'pushbutton';
        Formats(end,1).size = [140 25];
        Formats(end,1).callback = {@add_data,'','@selectmeasures','@selectsubjects'};
        Formats(end,1).span = [1 2];
    end
    
    [settings,Cancelled] = inputsdlg(Prompt,'Select Plot-input',Formats,settings,2);
    if Cancelled ~= 1
        if doavg == true
            if min(settings.selectsubjects) == 1
                doavg = true;
            else
                doavg = false;
            end
            settings.selectsubjects = setdiff(settings.selectsubjects-1,0);
            if isempty(settings.selectsubjects)
                settings.selectsubjects = 1:size(DATA,2);
            end
        end
        DATA_ALL = DATA;
        DATA = DATA(settings.selectmeasures,settings.selectsubjects);
        if doavg == true & size(DATA,2) > 1
            for i = 1:size(DATA,1)
                data = zeros(size(DATA,2),size(DATA(i,1).data,2));
                for j = 1:size(DATA,2)
                    data(j,:) = DATA(i,j).data;
                end
                data = mean(data,1);
                DATA(i,1).data = data;
                DATA(i,1).subject = 'Average';
                DATA(i,1).name = ['Average_' DATA(i,1).measure];
            end
            DATA = DATA(:,1);
        end
        if isfield(settings,'DATA_file')
            cfg.DATA_file = settings.DATA_file;
        end
    else
        DATA = [];
        DATA_ALL = [];
    end
    clearvars settings subjects measures Cancelled Prompt Formats
end
if ~isempty(DATA)
    cfg = lab_plot_match_numchans(cfg,1);
    cfg = lab_plot_check_numchans(cfg);
    cfg.PLOT = lab_plot_settings(cfg.PLOT,DATA,doIS);
end

% helper functions
    function add_graphmeasure(HandleM)
        if max(IsConnection) ~= 1
            return
        end
        connectivities = measures(IsConnection);
        
        Prompt2 = cell(0,2);
        Formats2 = [];
        
        Prompt2(end+1,:) = {'Select Connectivity data','selectconnectivity'};
        Formats2(end+1,1).type = 'list';
        Formats2(end,1).style = 'listbox';
        Formats2(end,1).format = 'input';
        Formats2(end,1).items = connectivities;
        Formats2(end,1).limits = [0 2]; % multi-select
        Formats2(end,1).size = [180 250];
        
        Prompt2(end+1,:) = {'Select Graph measure','selectgraph'};
        Formats2(end+1,1).type = 'list';
        Formats2(end,1).style = 'popupmenu';
        Formats2(end,1).format = 'input';
        Formats2(end,1).items = {'Degree','Betweenness centrality','Eigenvector centrality', ...
            'Clustering coefficient','Small dots','Large dots'};
        
        [settingsG,Cancelled2] = inputsdlg(Prompt2,'Add Graph measure',Formats2);
        if Cancelled2 == 1
            return
        end
        
        for SelectC = 1:length(settingsG.selectconnectivity)
            DATAnew = [];
            Nmeasure = find(strcmp(measures,settingsG.selectconnectivity{SelectC}),1);
            for SelectS = 1:size(DATA,2)
                matrix2 = lab_tril2matrix(DATA(Nmeasure,SelectS).data',ceil((2*length(DATA(Nmeasure,SelectS).data))^0.5));
                nodes = [];
                if strcmp(settingsG.selectgraph,'Degree')
                    if size(unique(matrix2),1) == 2 & min(matrix2) == 0 & max(matrix2) == 1
                        nodes = degrees_und(matrix2);
                    else
                        nodes = mean(matrix2,1);
                    end
                elseif strcmp(settingsG.selectgraph,'Betweenness centrality')
                    disp('    calculate betweenness centrality')
                    if size(unique(matrix2),1) == 2 & min(matrix2) == 0 & max(matrix2) == 1
                        [~,nodes] = edge_betweenness_bin(matrix2);
                        nodes = nodes ./ ((size(matrix2,1)-1)*(size(matrix2,1)-2));
                    else
                        matrixI = lab_invmatrix(matrix2);
                        [~,nodes] = edge_betweenness_wei(matrixI);
                        nodes = nodes ./ ((size(matrix2,1)-1)*(size(matrix2,1)-2));
                        clearvars matrixI
                    end
                elseif strcmp(settingsG.selectgraph,'Eigenvector centrality')
                    disp('    calculate eigenvector centrality')
                    nodes = eigenvector_centrality_und(matrix2);
                elseif strcmp(settingsG.selectgraph,'Clustering coefficient')
                    disp('    calculate Clustering coefficients')
                    if size(unique(matrix2),1) == 2 & min(matrix2) == 0 & max(matrix2) == 1
                        nodes = clustering_coef_bu(abs(matrix2))';
                    else
                        nodes = clustering_coef_wu(abs(matrix2))';
                    end
                elseif strcmp(settingsG.selectgraph,'Small dots')
                    nodes = ones(1,size(matrix2,1)) * 0.5;
                elseif strcmp(settingsG.selectgraph,'Large dots')
                    nodes = ones(1,size(matrix2,1));
                end
                if ~isempty(nodes)
                    nodes = nodes(:)';
                    DATAnew(1,SelectS).data = nodes; %#ok<AGROW>
                    DATAnew(1,SelectS).subject = DATA(Nmeasure,SelectS).subject; %#ok<AGROW>
                    DATAnew(1,SelectS).measure = [DATA(Nmeasure,SelectS).measure '_' regexprep(settingsG.selectgraph,' ','_')]; %#ok<AGROW>
                    DATAnew(1,SelectS).name = [DATAnew(1,SelectS).subject '_' DATAnew(1,SelectS).measure]; %#ok<AGROW>
                    DATAnew(1,SelectS).connections = false; %#ok<AGROW>
                    DATAnew(1,SelectS).nodes = true; %#ok<AGROW>
                    DATAnew(1,SelectS).labels = false; %#ok<AGROW>
                end
            end
            newmeasure = [DATA(Nmeasure,1).measure '_' regexprep(settingsG.selectgraph,' ','_')];
            if Nmeasure == length(measures)
                DATA = cat(1,DATA,DATAnew);
                measures = cat(1,measures,newmeasure);
                IsConnection = cat(1,IsConnection,false);
            else
                DATA = cat(1,DATA(1:Nmeasure,:),DATAnew,DATA(Nmeasure+1:end,:));
                measures = cat(1,measures(1:Nmeasure,:),newmeasure,measures(Nmeasure+1:end,:));
                IsConnection = cat(1,IsConnection(1:Nmeasure,:),false,IsConnection(Nmeasure+1:end,:));
            end
            set(HandleM,'String',measures,'Value',1:length(measures));
        end 
    end
    
    function add_data(HandleM,HandleS)
        disp('Matrix or File')
        Select2 = questdlg('Matrix or File?','Select Input','Matrix','File','File');
        DATAnew = lab_plot_read_datanew(settings,Labels,Select2);
        if isempty(DATAnew)
            return
        end
        if ~isempty(numchans)
            DATAnew = lab_plot_split_connections(DATAnew,[],numchans);
        end
        
        % select new measures
        Nm2 = strcmp(fieldnames(DATAnew),'measure');
        tmp2 = struct2cell(DATAnew);
        tmp2 = tmp2(Nm2,:,1);
        measuresnew = tmp2(:);
        if length(measuresnew) > 1
            Mselect = listdlg('PromptString','Measures:','SelectionMode','multiple', ...
                'ListString',measuresnew,'InitialValue',1:length(measuresnew),'CancelString', ...
                'None','ListSize',[200 400]);
            if isempty(Mselect)
                return
            end
            DATAnew = DATAnew(Mselect,:);
            measuresnew = measuresnew(Mselect);
        end
        clearvars tmp2 Mselect
        
        % match new and old data
        if size(DATAnew,2) > 1 & size(DATA,2) == 1
            Nm2 = strcmp(fieldnames(DATAnew),'subject');
            tmp2 = struct2cell(DATAnew);
            tmp2 = tmp2(Nm2,1,:);
            subjects2 = tmp2(:);
            Sselect = listdlg('PromptString','Files:','SelectionMode','single', ...
                'ListString',subjects2,'InitialValue',1,'CancelString', ...
                'None','ListSize',[200 400]);
            if isempty(Sselect)
                return
            end
            for m = 1:size(DATAnew,1)
                if strcmp(DATAnew(m,Sselect).measure,'Measure')
                    DATAnew(m,Sselect).measure = DATAnew(m,Sselect).subject;
                else
                    DATAnew(m,Sselect).measure = [DATAnew(m,Sselect).measure '_' DATAnew(m,Sselect).subject];
                end
                measuresnew{m} = DATAnew(m,Sselect).measure;
            end
            DATA = cat(1,DATA,DATAnew(:,Sselect));
        elseif size(DATA,2) > 1 | size(DATAnew,2) > 1
            Nm2 = strcmp(fieldnames(DATA),'subject');
            tmp2 = struct2cell(DATA);
            tmp2 = tmp2(Nm2,1,:);
            subjects1 = tmp2(:);
            Nm3 = strcmp(fieldnames(DATAnew),'subject');
            tmp3 = struct2cell(DATAnew);
            tmp3 = tmp3(Nm3,1,:);
            subjects2 = tmp3(:);
            [subjectsnew,Nsubj1,Nsubj2] = intersect(subjects1,subjects2);
            if isempty(Nsubj1) | isempty(Nsubj2)
                disp('Adding new data not possible, mismatch subjects/trials')
                return
            else
                DATA = cat(1,DATA(:,Nsubj1),DATAnew(:,Nsubj2));
                if ishandle(HandleS)
                    Smode = get(HandleS);
                    if strcmp(Smode.String{1},'--Average--')
                        subjectsnew = cat(1,cellstr('--Average--'),subjectsnew);
                        set(HandleS,'String',subjectsnew,'Value',1);
                    else
                        set(HandleS,'String',subjectsnew,'Value',1:length(subjectsnew));
                    end
                end
            end
        else
            DATA = cat(1,DATA,DATAnew);
        end
        measures = cat(1,measures,measuresnew);
        set(HandleM,'String',measures,'Value',1:length(measuresnew));
        IsConnection = find_connections(DATA);
    end
    
end

function IsConnection = find_connections(DATA)
    IsConnection = false(size(DATA,1),1);
    Nc = strcmp(fieldnames(DATA),'connections');
    if max(Nc) == 1
        tmp = struct2cell(DATA);
        tmp = tmp(Nc,:,1);
        tmp = tmp(:);
        for i = 1:length(tmp)
            if ~isempty(tmp{i}) & tmp{i} == true
                IsConnection(i) = true;
            end
        end
    end
end