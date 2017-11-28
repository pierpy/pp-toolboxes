% Helper script for lab_calculate_connectivity
%
% written by F.Hatz 2014

function cfg = lab_calculate_connect_writematrix(Result,cfg)

disp('    Write matrices')

if isempty(Result)
    return
end
if ~isfield(cfg,'Output_file') | ~isfield(cfg,'Output_filepath')
    return
end

Connect_filepath = cfg.Output_filepath;
if isfield(cfg,'patient') & ~isempty(cfg.patient)
    Connect_file = [cfg.patient '_Connectivity'];
else
    Connect_file = 'Connectivity';
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
    
    variables = fieldnames(result);
    tmp = [];
    for i = 1:size(variables,1)
        if isnumeric(result.(variables{i,1})) & size(result.(variables{i,1}),1) > 1 & ...
                size(result.(variables{i,1}),1) == size(result.(variables{i,1}),2) & ...
                ~strcmp(variables{i,1},'SL_hit') & ~strcmp(variables{i,1},'SLc_hit')
            tmp = [tmp i];
        end
    end
    if isempty(tmp)
        return
    else
        variables = variables(tmp,1);
    end
    
    for Nvar = 1:size(variables,1)
        if size(result.(variables{Nvar,1}),3) <= 100
            if isfield(result,'timestamp') & ~isempty(result.timestamp) & length(result.timestamp) == size(result.(variables{Nvar,1}),3)
                timestamp = result.timestamp;
            else
                for i = 1:size(result.(variables{Nvar,1}),3)
                    timestamp{i,1} = num2str(i);
                end
            end
            
            for nmatrix = 1:size(result.(variables{Nvar,1}),3)
                if ~isempty(Outputname)
                    matrixfile = fullfile(Connect_filepath,[Connect_file '_' Outputname '_' timestamp{nmatrix} '_' upper(variables{Nvar,1}) '_matrix.txt']);
                else
                    matrixfile = fullfile(Connect_filepath,[Connect_file(1:end-4) '_' timestamp{nmatrix} '_' upper(variables{Nvar,1}) '_matrix.txt']);
                end
                if ~exist(matrixfile,'file')
                    dlmwrite(matrixfile,result.(variables{Nvar,1})(:,:,nmatrix),'delimiter','\t','precision', 6);
                end
            end
        end
    end
end

end