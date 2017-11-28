% read data to TAPEEG format
%
% [data,header,cfg] = lab_read_data(Filename,cfg,nodata)
%
% cfg      structure with config (optional)
% nodata   true = only header
%
% Written by F.Hatz, Neurology Basel 2012

function [data,header,cfg] = lab_read_data(Filename,cfg,nodata,segment,novrb)

if ~exist('novrb','var')
    novrb = false;
end
if ~exist('segment','var')
    segment = []; %#ok<NASGU>
end
if ~exist('nodata','var')
    nodata = false; %#ok<NASGU>
end
if ~exist('Filename','var') | isempty(Filename);
    [Filename,Filepath]=uigetfile('*.*','Select file');
    if Filename ~= 0 & exist(fullfile(Filepath,'signal1.bin'),'file')
        disp('Found mff file')
        tmp=strfind(Filepath,filesep);
        Filename = [Filepath(tmp(end-1)+1:tmp(end) -1) '.mff'];
        Filepath = Filepath(1:tmp(end-1));
        clearvars tmp
    end
elseif exist(fullfile(Filename,'signal1.bin'),'file')
    disp('Found mff file')
    tmp=strfind(Filename,filesep);
    if tmp(end) == length(Filename)
        Filepath = Filename(1:tmp(end-1)-1);
        Filename = [Filename(tmp(end-1)+1:end-1) '.mff'];
    else
        Filepath = Filename(1:tmp(end)-1);
        Filename = [Filename(tmp(end)+1:end) '.mff'];
    end
    clearvars tmp
else
    [Filename,Filepath] = lab_filename(Filename);
    if exist(fullfile(Filepath,'signal1.bin'),'file')
        disp('Found mff file')
        tmp=strfind(Filepath,filesep);
        Filename = [Filepath(tmp(end-1)+1:end-1) '.mff'];
        Filepath = Filepath(1:tmp(end-1)-1);
        clearvars tmp
    end
end

if ~exist('cfg','var')
    cfg = [];
end

if isnumeric(Filename) & Filename == 0
    Filename = [];
    data = [];
    header = [];
elseif ~exist(fullfile(Filepath,Filename),'file') & ~exist([fullfile(Filepath,Filename(1:end-4)) filesep 'signal1.bin'],'file')
    Filename = [];
    data = [];
    header = [];
end

if ~isempty(Filename)
    [~,~,format] = lab_filename(Filename);
    format = lower(format);
    if exist(['lab_read_' format]) == 2 %#ok<EXIST>
        % script for reading files is selected by name (all scripts in folder 'iofile')
        eval(['tmp = nargout(@lab_read_' format ');']);
        if tmp > 2
            try
                eval(['tmp2 = nargin(@lab_read_' format ');']);
                if tmp2 >= 4
                    eval(['[data,header,cfg] = lab_read_' format '(fullfile(Filepath,Filename),cfg,nodata,segment);']);
                    header.sectionalreading = true;
                elseif tmp2 >= 3
                    eval(['[data,header,cfg] = lab_read_' format '(fullfile(Filepath,Filename),cfg,nodata);']);
                    header.sectionalreading = false;
                else
                    eval(['[data,header,cfg] = lab_read_' format '(fullfile(Filepath,Filename),cfg);']);
                    header.sectionalreading = false;
                end
                clearvars tmp2
                if isnumeric(data)
                    data = double(data); % TAPEEG needs double precision!
                end
            catch err
                if strcmp(format,'fif') | strcmp(format,'fiff')
                    try
                        data = lab_read_mri(fullfile(Filepath,Filename));
                        header = [];
                    catch err2
                        disp(getReport(err2))
                        data = [];
                        header = [];
                    end
                else
                    disp(getReport(err))
                    data = [];
                    header = [];
                end
            end
        else
            try
                eval(['[data] = lab_read_' format '(fullfile(Filepath,Filename));']);
                header = [];
                if ~isempty(data) & isnumeric(data) & size(data,1) == size(data,2)
                    header.datatype = 'matrix';
                    [patient,cfg] = lab_subjectname(fullfile(Matrix_filepath,Matrix_file),cfg);
                    header.subjects = cellstr(patient);
                end
            catch %#ok<CTCH>
                data = [];
                header = [];
            end
        end
        clearvars tmp
    else
        data =[];
        header = [];
        if isempty(Filename)
            disp('file not found')
        else
            % Try standard matlab importdata
            if ~isfield(cfg,'IMPORT') | ~isfield(cfg.IMPORT,'TXT') | isempty(cfg.IMPORT.TXT)
                cfg.IMPORT.TXT = false;
                Prompt = {'Import as txt-file','TXT'};
                Formats.type = 'check';
                [cfg.IMPORT,Cancelled] = inputsdlg(Prompt,'Import file',Formats,cfg.IMPORT);
                if Cancelled == 1 | isempty(cfg)
                    cfg.IMPORT.TXT = false;
                end
                clearvars Prompt Formats
            end
            if cfg.IMPORT.TXT == true
                [data,header,cfg] = lab_read_txt(fullfile(Filepath,Filename));
            else
                try
                    data = importdata(fullfile(Filepath,Filename));
                catch %#ok<CTCH>
                    disp('Importing file not possible')
                end
            end
        end
    end
end

if ~isempty(data) & ~isempty(header)
    if ~isfield(header,'filetype') | isempty(header.fileptye)
        if isstruct(data) & isfield(data,'anatomy')
            header.filetype = 'mri';
        elseif ~isfield(header,'filetype')
            if isnumeric(data) & ~isempty(data)
                if size(data,1) == size(data,2)
                    header.filetype = 'matrix';
                else
                    header.filetype = 'eeg';
                end
            end
        end
    end
    
    if isfield(cfg,'filename') & ~isempty(cfg.filename)
        [Filename,Filepath] = lab_filename(cfg.filename);
        cfg = rmfield(cfg,'filename');
    end
    cfg.EEG_file = Filename;
    cfg.EEG_filepath = Filepath;
    if ~isfield(header,'datatype') | ~strcmp(header.datatype,'matrix')
        [data,header,cfg] = lab_read_data_ext(data,header,cfg,novrb);
    end
elseif ~isempty(data) & ~isstruct(cfg)
    if ~isempty(Filename) & ~isnumeric(Filename)
        cfg = fullfile(Filepath,Filename);
    end
end

return;


