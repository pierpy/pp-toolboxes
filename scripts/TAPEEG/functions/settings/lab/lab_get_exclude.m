function exclude = lab_get_exclude(numchans)

global Main_Path ExcludeN Exclude
if exist(fullfile(Main_Path,'Exclude.xls'),'file') & isempty(Exclude)
    if ispc
        [ExcludeN,~,Exclude] = xlsread(fullfile(Main_Path,'Exclude.xls'),1);
    else
        [ExcludeN,~,Exclude] = xlsread(fullfile(Main_Path,'Exclude.xls'),1,'','basic');
    end
end

if ~isempty(find(ExcludeN==numchans,1)) & size(Exclude,2) == 2
    exclude = Exclude{find(ExcludeN==numchans,1),2};
    if ischar(exclude)
        exclude = str2num(exclude);
    end
else
    exclude = [];
end