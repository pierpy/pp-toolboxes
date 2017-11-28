function [bad,settings,skipautomated] = lab_get_filebad(header,settings,cfg)

if ~exist('header','var') | ~isfield(header,'EEG_file')
    header.EEG_file = [];
    header.EEG_filepath = [];
end
[~,~,~,filenameS_header] = lab_filename(header.EEG_file);
if ~exist('cfg','var') | ~isfield(cfg,'EEG_file')
    cfg.EEG_file = [];
    cfg.EEG_filepath = [];
end
[~,~,~,filenameS_cfg,filesegment] = lab_filename(cfg.EEG_file);
if isfield(header,'bad')
    bad = header.bad;
end
if ~exist('bad','var') | ~isfield(bad,'All')
    bad.All = [];
end

skipautomated = 0;
if ~strcmp(settings.filemethod,'skip')
    % Look for xls with bad channel information
    bad.file = [];
    bad.filename = '';
    if exist(fullfile(header.EEG_filepath,[filenameS_header '.xls']),'file')
        bad.file = eval_xls(fullfile(header.EEG_filepath,[filenameS_header '.xls']),filesegment);
        if ~isnan(bad.file)
            bad.filename = fullfile(header.EEG_filepath,[filenameS_header '.xls']);
            bad.filetxt = [filenameS_header '.xls'];
        else
            bad.file = [];
        end
    end
    if isempty(bad.filename) & exist(fullfile(header.EEG_filepath,[filenameS_header '.xlsx']),'file')
        bad.file = eval_xls(fullfile(header.EEG_filepath,[filenameS_header '.xlsx']),filesegment);
        if ~isnan(bad.file)
            bad.filename = fullfile(header.EEG_filepath,[filenameS_header '.xlsx']);
            bad.filetxt = [filenameS_header '.xlsx'];
        else
            bad.file = [];
        end
    end
    if isempty(bad.filename) & isfield(cfg,'EEG_filepath') & exist(fullfile(cfg.EEG_filepath,[filenameS_cfg '.txt']),'file')
        bad.file = eval_txt(fullfile(cfg.EEG_filepath,[filenameS_cfg '.txt']));
        if ~isnan(bad.file)
            bad.filename = fullfile(cfg.EEG_filepath,[filenameS_cfg '.txt']);
            bad.filetxt = [filenameS_cfg '.txt'];
        else
            bad.file = [];
        end
    end
    if isempty(bad.filename) & isfield(cfg,'EEG_filepath') & exist(fullfile(cfg.EEG_filepath,[filenameS_cfg '.txt']),'file')
        bad.file = eval_txt(fullfile(cfg.EEG_filepath,[filenameS_header '.txt']));
        if ~isnan(bad.file)
            bad.filename = fullfile(cfg.EEG_filepath,[filenameS_header '.txt']);
            bad.filetxt = [filenameS_header '.txt'];
        else
            bad.file = [];
        end
    end
    if isempty(bad.filename) & isfield(header,'badchans')
        bad.file = header.badchans;
        bad.filename = fullfile(cfg.EEG_filepath,[filenameS_header '.info']);
         bad.filetxt = [filenameS_header '.info'];
    end
    if ~isempty(bad.filename)
        bad.All = union(bad.All,bad.file);
        if strcmp(settings.filemethod,'or')
            skipautomated = 1;
        end
    else
        bad.filetxt = 'Missing file or wrong file-format';
    end
    if strcmp(settings.filemethod,'only')
        skipautomated = 1;
    end
end

end

function bad = eval_xls(Filename,filesegment)

if ispc
    [~,~,badxls] = xlsread(Filename);
else
    [~,~,badxls] = xlsread(Filename,1,'','basic');
end
if length(badxls{1,1}) >= 12 & (strcmp(badxls{1,1}(1:12),'Bad channels') | strcmp(badxls{1,1}(1:12),'bad channels'))
    disp('   Select bad channels by xls-file')
    filenr = filesegment;
    if filenr > (size(badxls,1)-1)
        filenr = 2;
    elseif isempty(filenr)
        filenr = 2;
    else
        filenr = filenr + 1;
    end
    tmp = badxls{filenr,1};
    if ischar(tmp)
        tmp = str2num(tmp); %#ok<ST2NM>
    end
    if length(tmp) == 1 & (tmp == 0 | isnan(tmp))
        bad = [];
    else
        bad = sort(tmp);
    end
else
    disp('    xls-file has no bad channels information (check header)')
    bad = NaN;
end

end

function bad = eval_txt(Filename)

fid=fopen(Filename,'r');
tline=fgetl(fid);
if strcmp(tline,'Excluded channels:')
    disp('   Select bad channels by txt-file')
    bad = str2num(fgetl(fid)); %#ok<ST2NM>
    if length(bad) == 1 & (bad == 0 | isnan(bad))
        bad = [];
    else
        bad = sort(bad);
    end
else
    disp('    txt-file has no bad channels information')
    bad = NaN;
end
fclose(fid);

end