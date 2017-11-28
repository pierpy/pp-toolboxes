% result = lab_read_csv(filename)
% 
% by F. Hatz, Neurology Basel

function result = lab_read_csv(Filename)

if ~exist('Filename','var')
    [Filename,Filepath] = uigetfile('*.csv','Select csv-file');
    Filename = fullfile(Filepath,Filename);
end

try
    result = importdata(Filename);
catch %#ok<CTCH>
    result = [];
    return
end

if ~isempty(result)
    [~,Filepath,~,FilenameS] = lab_filename(Filename);
    if exist(fullfile(Filepath,[FilenameS '.head']),'file')
        try %#ok<TRYNC>
            load(fullfile(Filepath,[FilenameS '.head']),'-mat');
        end
        if exist('header','var') & isstruct(header) %#ok<NODEF>
            if isfield(header,'subjects') & ~isempty(header.subjects) & length(header.subjects) == size(result,2) & ...
                    isfield(header,'measures') & ~isempty(header.measures) & length(header.measures) == size(result,1)
                result = cat(1,header.subjects(:)',num2cell(result));
                if ~isfield(header,'main') | isempty(header.main)
                    header.main = {''};
                end
                result = cat(2,cat(1,header.main(1),header.measures(:)),result);
            elseif isfield(header,'subjects') & ~isempty(header.subjects) & length(header.subjects) == size(result,2)
                result = cat(1,header.subjects(:)',num2cell(result));
            elseif isfield(header,'measures') & ~isempty(header.measures) & length(header.measures) == size(result,1)
                result = cat(2,header.measures(:),num2cell(result));
            end
            return
        end
    end    
    if length(FilenameS) > 9 & strcmp(FilenameS(end-9:end),'_bootstrap') & ...
            exist(fullfile(Filepath,[FilenameS(1:end-10) '.xlsx']),'file')
        Filename2 = fullfile(Filepath,[FilenameS(1:end-10) '.xlsx']);
    elseif length(FilenameS) > 9 & strcmp(FilenameS(end-9:end),'_bootstrap') & ...
            exist(fullfile(Filepath,[FilenameS(1:end-10) '.xls']),'file')
        Filename2 = fullfile(Filepath,[FilenameS(1:end-10) '.xls']);
    else
        return
    end
    try
        [~,~,header] = xlsread(Filename2,1);
    catch %#ok<CTCH>
        return
    end
    if isempty(header)
        return
    end
    if size(header,1) == size(result,1)+1 & size(header,2) == size(result,2)+1
        result = cat(2,header(2:end,1),num2cell(result));
        result = cat(1,header(:,1),result);
    elseif size(header,1) == size(result,2)+1 & size(header,2) == size(result,1)+1
        result = cat(2,header(1,2:end)',num2cell(result));
        result = cat(1,header(:,1)',result);
        result = result';
    elseif size(header,1) == size(result,1)
        result = cat(2,header(:,1),num2cell(result));
    elseif size(header,2) == size(result,2)
        result = cat(1,header(1,:),num2cell(result));
    elseif size(header,1) == size(result,2)
        result = cat(1,header(:,1)',num2cell(result));
        result = result';
    elseif size(header,2) == size(result,1)
        result = cat(2,header(1,:)',num2cell(result));
        result = result';
    elseif size(header,1) == size(result,1)+1
        result = cat(2,header(2:end,1),num2cell(result));
    elseif size(header,2) == size(result,2)+1
        result = cat(1,header(1,2:end),num2cell(result));
    elseif size(header,1) == size(result,2)+1
        result = cat(1,header(2:end,1)',num2cell(result));
        result = result';
    elseif size(header,2) == size(result,1)+1
        result = cat(2,header(1,2:end)',num2cell(result));
        result = result';
    end
end


