% For result-files with many variables per measure (e.g. EEG-channels),
% this function will use a Mappings-file (*.xls) to group variables in
% mappings and output the reult as a new Excel-file (_xls2map.xls)
%
% lab_xls2mappings(xlsinput,mappings)
%
% written by F. Hatz 2014

function lab_xls2mappings(xlsinput,mappings)

% Read xls-file
if ~exist('xlsinput','var') | isempty(xlsinput)
    [xlsinput,~,~,~,cfg] = lab_read_statistics([],-1,0,1,0,0);
end

if ~isempty(xlsinput)
    header.file = cfg.file;
    header.path = cfg.path;
    header.subjects = xlsinput(1,2:end);
    header.vars = xlsinput(2:end,1);
    try
        tmp = xlsinput(2:end,2:end);
        tmp2 = find(cellfun(@isempty,tmp));
        if ~isempty(tmp2)
            for i = 1:length(tmp2);
                tmp{tmp2(i)} = NaN;
            end
        end
        data = cell2mat(tmp);
    catch %#ok<CTCH>
        disp('Input file has invalid / non numeric data')
        return
    end
    clearvars xlsinput tmp tmp2
else
    return
end
if isempty(data)
    return
end

if ~isfield(cfg,'numresults')
    cfg.numresults = 0;
end

if ~exist('mappings','var')
    cfg.Mappings = [];
else
    cfg.Mappings = mappings;
end

Prompt = {};
Formats = {};

Prompt(end+1,:) = {'Number of variables per measure', 'clustervars'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'integer';
Formats(end,1).limits = [0 9999];
Formats(end,1).size = 60;

Prompt(end+1,:) = {'Mappings', 'Mappings'};
Formats(end+1,1).type = 'check';
Formats(end,1).callback = {@lab_load_mappings,'Mappings','Mappings'};

[cfg,Cancelled] = inputsdlg(Prompt,'Excel to Mappings',Formats,cfg);
if Cancelled == 1
    return
end

if isempty(cfg.Mappings) | ~isfield(cfg.Mappings,'mappingsChannels') | cfg.Mappings.mappingsChannels ~= cfg.clustervars
    disp('Abort: Mappings-file and Inputdata not matching')
    return
end

if cfg.clustervars < 1 | mod(size(data,1)-cfg.numresults,cfg.clustervars) > 0
    return
else
    cfg.numclusters = (size(data,1)-cfg.numresults) / cfg.clustervars;
end

if cfg.numresults > 0
    results = data(end-cfg.numresults+1:end,:);
    data = data(1:end-cfg.numresults,:);
    header.results = header.vars(end-cfg.numresults+1:end,1);
    header.vars = header.vars(1:end-cfg.numresults,1);
else
    results = [];
    header.results = {};
end

varnames = cell(0,1);
for i = 1:cfg.numclusters
    tmp = strfind(header.vars{(i-1)*cfg.clustervars+1,1},'_');
    if ~isempty(tmp)
        measure = header.vars{(i-1)*cfg.clustervars+1,1}(1:tmp(1)-1);
    else
        measure = ['Measure' num2str(i)];
    end
    measure = regexprep(measure,{'::',':',';','/','\',','},'_');
    for j = 1:size(cfg.Mappings.mappingstitle,1)
        varnames{end+1,1} = [measure '_' cfg.Mappings.mappingstitle{j,1}];
    end
    clearvars measure
end
clearvars tmp i j

dataout = [];
datatmp = zeros(size(cfg.Mappings.mappings,2),size(data,2));
for i = 1:cfg.numclusters
    dataC = data((i-1)*cfg.clustervars+1:i*cfg.clustervars,:);
    for j = 1:size(cfg.Mappings.mappings,2)
        datatmp(j,:) = mean(dataC(cfg.Mappings.mappings{1,j},:),1);
    end
    dataout = cat(1,dataout,datatmp);
end
dataout = cat(1,dataout,results);
varnames = cat(1,varnames,header.results);

xlsout = [cellstr(['C' num2str(size(cfg.Mappings.mappings,2)) ' R' num2str(cfg.numresults)]) header.subjects];
xlsout = cat(1,xlsout,cat(2,varnames,num2cell(dataout)));

if isfield(header,'path')
    if size(xlsout,2) > 255
        fileout = fullfile(header.path,[header.file '_xls2map.xlsx']);
    else
        fileout = fullfile(header.path,[header.file '_xls2map.xls']);
    end
    lab_write_xls(fileout,xlsout);
    clearvars xlsout tmp
end

return