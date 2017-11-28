function datainput = lab_set_NaN(datainput)
    
tmp = ~cellfun(@isnumeric,datainput);
tmp2 = cellfun(@isempty,datainput);
tmp(tmp2) = true;
if sum(tmp(:,1)) == size(tmp,1);
    flag = true;
else
    flag = false;
end
tmp(sum(tmp,2) == size(tmp,2),:) = false;
if flag == true
    tmp(:,1) = false;
end
datainput(tmp) = {NaN};

return