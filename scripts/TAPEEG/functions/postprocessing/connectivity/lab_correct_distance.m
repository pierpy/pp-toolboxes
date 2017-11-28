% Helper script for lab_calculate_connectivity
%
% written by F.Hatz 2014

function matrix = lab_correct_distance(matrix,distance)

if isnumeric(matrix)
    for nmatrix = 1:size(matrix,3)
        Mtmp = matrix(:,:,nmatrix);
        tmp = Mtmp;
        tmp(1:size(tmp,1)+1:end) = 0;
        tmp = max(tmp(:));
        Mtmp(distance==1) = tmp;
        matrix(:,:,nmatrix) = Mtmp;
    end
elseif isstruct(matrix)
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
    clearvars tmp
    for Nvar = 1:size(variables,1)
        disp(['    correct for nearest neighboors ' variables{Nvar,1}])
        for nmatrix = 1:size(matrix.(variables{Nvar,1}),3)
            Mtmp = matrix.(variables{Nvar,1})(:,:,nmatrix);
            tmp = Mtmp;
            tmp(1:size(tmp,1)+1:end) = 0;
            tmp = max(tmp(:));
            Mtmp(distance==1) = tmp;
            matrix.(variables{Nvar,1})(:,:,nmatrix) = Mtmp;
        end
    end
end