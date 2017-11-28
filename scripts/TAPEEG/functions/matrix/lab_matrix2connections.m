function [xlsout,cfg] = lab_matrix2connections(cfg,flagstore)

if ~exist('flagstore','var')
    flagstore = true;
end
if ~exist('cfg','var')
    cfg = [];
end
if isstruct(cfg) & isfield(cfg,'data_file')
    [matrix,header,cfg] = lab_read_matrix(fullfile(cfg.data_filepath,cfg.data_file));
elseif ~isempty(cfg) & isnumeric(cfg) & size(cfg,1) == size(cfg,2)
    matrix = cfg;
    header = [];
    cfg = [];
else
    cfg.SEARCH.searchstring = {'matrix.txt'};
    [Filelist,cfg] = lab_search_files(cfg,0,1);
    Filelist = Filelist.Filelist;
    if isempty(Filelist)
        return
    end
    [matrix,header,cfg] = lab_read_matrix(Filelist,cfg);
    cfg.EEG_filepath = cfg.SEARCH.searchfolder;
end
if size(matrix,1) ~= size(matrix,2)
    xlsout = [];
    disp('No valid input-file for matrix')
    return
end
if isfield(header,'channels')
    labels = cellstr(header.channels);
else
    labels = cell(size(matrix,1),1);
    for i = 1:size(matrix,1)
        labels{i,1} = ['Ch' num2str(i,'%03d')];
    end
end
[connections,labels] = lab_extract_tril(matrix,labels,1);
if isfield(header,'subjects') & length(header.subjects) == size(connections,2)
    subjects = header.subjects;
else
    subjects = cell(1,size(connections,2));
    for i = 1:size(connections,2)
        subjects{1,i} = ['Trial' num2str(i,'%02d')];
    end
end
xlsout = cat(2,{['C' num2str(size(connections,1)) ' R0']},subjects);
xlsout = cat(1,xlsout,cat(2,labels,num2cell(connections)));
if isfield(cfg,'EEG_filepath') & exist(cfg.EEG_filepath,'dir') & flagstore == true
    [~,~,~,Filename] = lab_filename(cfg.EEG_file);
    lab_write_xls(fullfile(cfg.EEG_filepath,[Filename '.xlsx']),xlsout);
end

end