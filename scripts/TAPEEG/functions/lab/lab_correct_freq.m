function string = lab_correct_freq(string)

[~,~,~,string] = lab_filename(string);
tmp = regexpi(string,'_F.*_.*_');
if isempty(tmp)
    return
end
tmp2 = strfind(string(tmp(1)+1:end),'_');
if isempty(tmp2)
    return
end
if length(string) < tmp(1)+tmp2(1)+2 | isempty(str2double(string(tmp(1)+2:tmp(1)+tmp2(1)-1))) | ...
        isempty(str2double(string(tmp(1)+tmp2(1)+1)))
    return
end
string(tmp(1)+tmp2(1)) = 'F';