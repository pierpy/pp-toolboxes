% Function to perform kmeans-analysis for results of connectivity analysis
%
% Result = lab_calculate_connect_kmeans(Result,cfg)
%
% Written by F. Hatz 2014 Neurology Basel

function Result = lab_calculate_connect_kmeans(Result,cfg)

disp('Calculate Kmeans-Clustering')

if isempty(Result) | ~isfield(cfg,'CONNECT') | ~isfield(cfg.CONNECT,'clustering') | ...
        isempty(cfg.CONNECT.clustering) | cfg.CONNECT.clustering <= 1
    return
end

frequencies = fieldnames(Result);
tmp = [];
for i = 1:size(frequencies,1)
    if isstruct(Result.(frequencies{i,1})) & strcmp(frequencies{i,1}(1),'F')
        tmp = [tmp i];
    end
end
if isempty(tmp)
    frequencies = 0;
else
    frequencies = frequencies(tmp,1);
end

for Nfreq = 1:size(frequencies,1)
    if ~isnumeric(frequencies)
        result = Result.(frequencies{Nfreq,1});
        Outputname = frequencies{Nfreq,1};
    else
        result = Result;
        Outputname = '';
    end
    
    if isfield(cfg,'Output_filepath')
        Connect_filepath = cfg.Output_filepath;
        if isfield(cfg,'patient') & ~isempty(cfg.patient)
            Connect_file = [cfg.patient '_Kmeans'];
            patient = cfg.patient;
        else
            Connect_file = 'Kmeans';
            patient = [];
        end
        if ~isempty(Outputname)
            Connect_file = [Connect_file '_' Outputname '.mat'];
        else
            Connect_file = [Connect_file '.mat'];
        end
    end
    
    variables = fieldnames(result);
    tmp = [];
    for i = 1:size(variables,1)
        if isnumeric(result.(variables{i,1})) & size(result.(variables{i,1}),1) > 1 & ...
                size(result.(variables{i,1}),1) == size(result.(variables{i,1}),2) & ...
                ~strcmp(variables{i,1},'SL_hit') & ~strcmp(variables{i,1},'SLc_hit') & ...
                ~strcmp(variables{i,1},'wplv_angle')
            tmp = [tmp i];
        end
    end
    if isempty(tmp)
        return
    else
        variables = variables(tmp,1);
    end
    
    for Nvar = 1:size(variables,1)
        tmp = lab_calculate_kmeans(result.(variables{Nvar,1}),cfg.CONNECT.clustering);
        result.(['kmeans_' variables{Nvar,1}]) = tmp;
    end
    
    if exist('Connect_file','var')
        save(fullfile(Connect_filepath,Connect_file),'result','patient');
    end
    
    if ~isnumeric(frequencies)
        Result.(frequencies{Nfreq,1}) = result;
    else
        Result = result;
    end
end

end