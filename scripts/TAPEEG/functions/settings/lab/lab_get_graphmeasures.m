function settings = lab_get_graphmeasures(settings)

if ~exist('settings','var')
    settings = [];
end

if max(strcmp(settings.measure,'ClusteringCoeff_Norm')) | max(strcmp(settings.measure,'ShortestPath_Norm'))
    clearvars Prompt Formats
    dodiag = 0;
    if ~isfield(settings,'randnumber') | isempty(settings.randnumber)
        settings.randnumber = 50;
        settings.randiter = 5;
        dodiag = 1;
    end
    
    if dodiag == 1
        Prompt = cell(0,2);
        Formats = [];
        
        Prompt(end+1,:) = {'Number of random matrices','randnumber'};
        Formats(end+1,1).type = 'edit';
        Formats(end,1).format = 'integer';
        Formats(end,1).limits = [0 999];
        Formats(end,1).size = 30;
        
        Prompt(end+1,:) = {'Iterations for random matrices','randiter'};
        Formats(end+1,1).type = 'edit';
        Formats(end,1).format = 'integer';
        Formats(end,1).limits = [0 999];
        Formats(end,1).size = 30;
        
        [settings,Cancelled] = inputsdlg(Prompt,'Random matrices',Formats,settings);
        if Cancelled == 1
            settings.randnumber = [];
            tmp = ~strcmp(settings.measure,'ClusteringCoeff_Norm');
            settings.measure = settings.measure(tmp);
            tmp = ~strcmp(settings.measure,'ShortestPath_Norm');
            settings.measure = settings.measure(tmp);
        end
    end
else
    settings.randnumber = [];
    settings.randiter = [];
    if ~isempty(settings.measure)
        tmp = ~strcmp(settings.measure,'ClusteringCoeff_Norm');
        settings.measure = settings.measure(tmp);
        tmp = ~strcmp(settings.measure,'ShortestPath_Norm');
        settings.measure = settings.measure(tmp);
    end
end

if max(strcmp(settings.measure,'MST-Ref'))
    if ~isfield(settings,'MSTref') | isempty(settings.MSTref)
        settings.MSTref = lab_load_matrix(settings.MSTref,[],true);
    end
    if isempty(settings.MSTref)
        tmp = ~strcmp(settings.measure,'MST-Ref');
        settings.measure = settings.measure(tmp);
    end
else
    settings.MSTref = [];
    tmp = ~strcmp(settings.measure,'MST-Ref');
    if ~isempty(settings.measure)
        settings.measure = settings.measure(tmp);
    end
end

if max(strcmp(settings.measure,'Matrix-Ref'))
    if ~isfield(settings,'MatrixRef') | isempty(settings.MatrixRef)
        settings.MatrixRef = lab_load_matrix(settings.MatrixRef);
    end
    if isempty(settings.MatrixRef)
        tmp = ~strcmp(settings.measure,'Matrix-Ref');
        settings.measure = settings.measure(tmp);
    end
else
    settings.MatrixRef = [];
    tmp = ~strcmp(settings.measure,'Matrix-Ref');
    if ~isempty(settings.measure)
        settings.measure = settings.measure(tmp);
    end
end

pause(0.2);