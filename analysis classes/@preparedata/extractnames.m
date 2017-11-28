function [ obj ] = extractnames ( obj, path )
    %EXTRACT_NAMES returns a cell whit names of files in path (specificare full path)
    f = dir(path);
    f(1)=[];
    f(1)=[];
    isD = [f.isdir];
    f = f(~isD);
    obj.filesName={};
    for i = 1:size(f,1)
        obj.filesName{i}=f(i).name;
    end
end