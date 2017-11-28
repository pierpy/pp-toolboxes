% Helper file for lab_plot_IS
%
% Written by F. Hatz 2013

function [matrixedges,result] = lab_calc_edges(numchans)

result = 1:numchans*numchans;
result = reshape(result,numchans,numchans);
result(1:numchans+1:end) = 0;

for i = 1:numchans*numchans;
    [y{i,1},x{i,1}] = find(result == i); %#ok<AGROW>
end
value = num2cell(1:numchans*numchans)';
result2 = [x y value];
clearvars x y value

tmp = ~cellfun('isempty',result2(:,1));
result2 = result2(tmp,:);
clearvars tmp

list = [];
j = 1;
for i = 1:(numchans*(numchans-1))
    tmp = cell2mat(result2(i,1:2));
    tmp = num2str(sort(tmp));
    tmp2 = find(strcmp(list,tmp));
    if isempty(tmp2)
        tmp = cellstr(tmp);
        list = [list,tmp]; %#ok<AGROW>
        matrixedges(j,1:3) = cell2mat(result2(i,1:3)); %#ok<AGROW>
        j = j+1;
    else
        matrixedges(tmp2,4) = result2{i,3}; %#ok<AGROW>
    end
    clearvars tmp tmp2
end
clearvars list i j result2

for i = 1:(numchans*(numchans-1))/2
    result(matrixedges(i,1),matrixedges(i,2)) = matrixedges(i,3);
    result(matrixedges(i,2),matrixedges(i,1)) = matrixedges(i,4);
end
        
    