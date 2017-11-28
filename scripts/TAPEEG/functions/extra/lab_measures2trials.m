function lab_measures2trials

[data,header,~,~,cfg] = lab_read_statistics([],0,0,1,0,1);

if isfield(header,'variables') & mod(size(data,2),cfg.clustervars) == 0
    vars = header.variables;
    measures = header.measures;
else
    vars = header.measures;
    measures = {'Measure'};
    cfg.clustervars = length(vars);
    cfg.numclusters = 1;
end
subjects = header.subjects;

xlsout = {''};
for i = 1:length(subjects)
    for j = 1:length(vars)
        if length(subjects) == 1
            xlsout{1,end+1} = vars{j}; %#ok<AGROW>
        else
            xlsout{1,end+1} = [subjects{i} '_' vars{j}]; %#ok<AGROW>
        end
    end
end

group = {'group'};
for i = 1:size(data,1);
    group = [group num2cell(1:cfg.clustervars)]; %#ok<AGROW>
end
data = reshape(data,[size(data,1) cfg.clustervars cfg.numclusters]);
data = permute(data,[2 1 3]);
data = reshape(data,size(data,1)*size(data,2),size(data,3));

xlsout = cat(1,xlsout,[measures num2cell(data')]);
xlsout = cat(1,xlsout,group);
xlsout{1,1} = 'C1 R0';

if strcmp(cfg.file(end-10:end),'_statistics')
    cfg.file = cfg.file(1:end-11);
end

[Filename,Filepath] = uiputfile([cfg.file '_M2T.xlsx'],'Select File to store result');
if ~ischar(Filename) | Filename == 0
    return
end
lab_write_xls(fullfile(Filepath,Filename),xlsout);

end
        