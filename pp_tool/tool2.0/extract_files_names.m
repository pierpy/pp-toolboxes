function [ names_cell ] = extract_files_names( path )
%EXTRACT_NAMES returns a cell whit names of files in path (specificare full path)
    f = dir(path);
    f(1)=[];
    f(1)=[];
    isD = [f.isdir];
    f = f(~isD);
    names_cell={};
    for i = 1:size(f,1)
        names_cell{i,1}=f(i).name;
    end

end
