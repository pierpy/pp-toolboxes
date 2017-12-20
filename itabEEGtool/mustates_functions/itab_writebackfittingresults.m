function itab_writebackfittingresults(filename, VariablesToWrite, subjects )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
    vars = fieldnames(VariablesToWrite);
    microstates = {'A', 'B', 'C', 'D'};
    T = [];
    for i = 1: length(vars) 
        eval(strcat('varr = ','VariablesToWrite.', vars{i},'.'''))
        T = [T, varr];
    end
    VarNames = {};
    offset = 0;
    for j = 1: length(vars)
        if ~strcmp(vars{j}(1:end-2), 'TotalGev')
        for k = 1: 4
            offset = offset + 1;
            VarNames{offset} = strcat(vars{j}, '_', microstates{k});
        end
        else
            offset = offset + 1;
            VarNames{offset} = vars{j};
        end
    end
    data=T;     %Sample 2-dimensional data
    data_cells=num2cell(data);     %Convert data to cell array
    col_header=VarNames;     %Row cell array (for column labels)
    row_header=subjects;     %Column cell array (for row labels)
    output_matrix=[{' '} col_header; row_header data_cells];     %Join cell arrays
    xlswrite(filename,output_matrix);     %Write data and both headers
end

