% Function to perform graph analysis for results of connectivity analysis
%
% [Result,cfg] = lab_calculate_connect_dograph(Result,header,cfg)
%
% Written by F. Hatz 2014 Neurology Basel

function [Result,cfg] = lab_calculate_connect_dograph(Result,header,cfg)

if isempty(Result)
    return
end

if isfield(cfg,'Output_file')
    Output_file = cfg.Output_file;
    Output_fileS = cfg.Output_fileS;
    Output_filepath = cfg.Output_filepath;
end
if isfield(cfg,'patient')
    Graph_file = [cfg.patient '_'];
else
    Graph_file = '';
end

frequencies = fieldnames(Result);
tmp = [];
for i = 1:size(frequencies,1)
    if isstruct(Result.(frequencies{i,1})) & strcmp(frequencies{i,1}(1),'F')
        tmp = [tmp i]; %#ok<AGROW>
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
    
    variables = fieldnames(result);
    tmp = [];
    for i = 1:size(variables,1)
        if isnumeric(result.(variables{i,1})) & size(result.(variables{i,1}),1) > 1 & ...
                size(result.(variables{i,1}),1) == size(result.(variables{i,1}),2) & ...
                ~strcmp(variables{i,1},'SL_hit') & ~strcmp(variables{i,1},'SLc_hit') & ...
                ~strcmp(variables{i,1},'wplv_angle')
            tmp = [tmp i]; %#ok<AGROW>
        end
    end
    if isempty(tmp)
        return
    else
        variables = variables(tmp,1);
    end
    
    for Nvar = 1:size(variables,1)
        if strcmp(variables{Nvar,1},'dpli')
            cfg.isdpli = 1;
        else
            cfg.isdpli = 0;
        end
        if ~isempty(Outputname)
            cfg.Output_file = [Graph_file upper(variables{Nvar,1}) '_' Outputname '.xls'];
        else
            cfg.Output_file = [Graph_file upper(variables{Nvar,1}) '.xls'];
        end
        if isfield(result,[variables{Nvar,1} '_timestamp']) & ~isempty(result.([variables{Nvar,1} '_timestamp'])) & ...
                length(result.([variables{Nvar,1} '_timestamp'])) == size(result.(variables{Nvar,1}),3)
            [result.([variables{Nvar,1} '_graph']),cfg.CONNECT,cfg] = lab_graphanalysis(result.(variables{Nvar,1}),header,cfg.CONNECT,cfg,result.([variables{Nvar,1} '_timestamp']),1);
        else
            [result.([variables{Nvar,1} '_graph']),cfg.CONNECT,cfg] = lab_graphanalysis(result.(variables{Nvar,1}),header,cfg.CONNECT,cfg,[],1);
        end
        cfg = rmfield(cfg,'isdpli');
    end
    
    if ~isnumeric(frequencies)
        Result.(frequencies{Nfreq,1}) = result;
    else
        Result = result;
    end
end

if exist('Output_file','var')
    cfg.Output_file = Output_file;
    cfg.Output_fileS = Output_fileS;
    cfg.Output_filepath = Output_filepath;
elseif isfield(cfg,'Output_file')
    cfg = rmfield(cfg,'Output_file');
end

end