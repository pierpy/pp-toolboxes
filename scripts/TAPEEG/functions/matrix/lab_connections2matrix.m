function Result = lab_connections2matrix

Result = [];

[data_file,data_filepath]=uigetfile('*.xls;*.xlsx','Select xls-file with connections');
if ~ischar(data_file) | isempty(data_file)
    return
end
datainput = lab_read_xls(fullfile(data_filepath,data_file));
datainput = lab_correctheader(datainput);
[datainput,cfg] = lab_getstructure(datainput);

Vars = datainput(2:end,1);
Subjects = datainput(1,2:end);
try
    connections = cell2mat(datainput(2:end,2:end));
catch %#ok<CTCH>
    return
end
if cfg.clustervars > 1
    Nmeasures = floor(size(connections,1) / cfg.clustervars);
    Vars = Vars(1:Nmeasures*cfg.clustervars,:);
    connections = connections(1:Nmeasures*cfg.clustervars,:);
    button = questdlg('Diagonal included?','Diagonal included?','No','Yes','Yes');
    if strcmp(button,'Yes')
        dodiag = true;
        if round(cfg.clustervars^0.5) == cfg.clustervars^0.5
            Nchans = cfg.clustervars^0.5;
            matrix = reshape(connections,[Nchans Nchans Nmeasures*size(connections,2)]);
        else
            Nchans = floor((2 * cfg.clustervars)^0.5);
            matrix = lab_tril2matrix(connections,Nchans);
        end
    else
        dodiag = false;
        Nchans = ceil((2 * cfg.clustervars)^0.5);
        matrix = lab_tril2matrix(connections,Nchans);
    end
else
    Nchans = inputdlg('Number of channels');
    if isempty(Nchans)
        return
    end
    Nchans = str2num(Nchans{1}); %#ok<ST2NM>
    if mod(size(connections,1),Nchans^2/2 + Nchans/2) == 0
        dodiag = true;
    else
        dodiag = false;
    end
    matrix = lab_tril2matrix(connections,Nchans);
end

if isempty(matrix)
    disp('Abort: wrong input data')
    return
end

Labels = cell(size(matrix,1)+1,1);
Ltmp = cell(length(Vars),1);
Mtmp = cell(length(Vars),1);
for i = 1:length(Vars);
    tmp = strfind(Vars{i},'_');
    if ~isempty(tmp)
        Ltmp{i} = Vars{i}(tmp(end)+1:end);
        Mtmp{i} = Vars{i}(1:tmp(end)-1);
    else
        Ltmp{i} = Vars{i};
        Mtmp{i} = Vars{i};
    end
end
for i = 1:length(Labels)-1
    tmp = strfind(Ltmp{i},'-');
    if ~isempty(tmp)
        if i == 1
            Labels{i} = Ltmp{i}(1:tmp(1)-1);
        end
        Labels{i+1} = Ltmp{i}(tmp(1)+1:end);
    else
        tmp = floor(length(Ltmp{i}) / 2);
        if i == 1
            Labels{i} = Ltmp{i}(1:tmp);
        end
        Labels{i+1} = Ltmp{i}(tmp+1:end);
    end
end
if dodiag == true
    Labels = Labels(2:end);
else
    Labels = Labels(1:end-1);
end

Measures = unique(Mtmp,'stable');
NumM = length(Measures);
NumS = length(Subjects);
if NumS * NumM ~= size(matrix,3);
    NumM = size(matrix,3) / NumS;
    if NumM ~= round(NumM)
        disp('Abort: wrong input data')
        return
    end
    Measures = cell(NumM,1);
    for i = 1:NumM
        Measures{i} = ['M' num2str(i)];
    end
end
OutputName = {};
for i = 1:NumS
    for j = 1:NumM
        OutputName{1,end+1} = [Subjects{i} '_' Measures{j}]; %#ok<AGROW>
    end
end
Output_path = fullfile(data_filepath,'Matrices');
warning off %#ok<WNOFF>
mkdir(Output_path);
warning on %#ok<WNON>
for i = 1:length(OutputName)
    lab_write_matrix(fullfile(Output_path,[OutputName{i} '_matrix.txt']),matrix(:,:,i));
end

matrix = reshape(matrix,[size(matrix,1) size(matrix,2) NumM NumS]);
matrix = permute(matrix,[1 2 4 3]);
for i = 1:NumM
    Result.(Measures{i}).matrix = matrix(:,:,:,i);
    Result.(Measures{i}).subjects = Subjects;
    Result.(Measures{i}).channels = char(Labels);
end
[~,~,~,data_fileS] = lab_filename(data_file);
save(fullfile(Output_path,[data_fileS '.mat']),'Result');

end