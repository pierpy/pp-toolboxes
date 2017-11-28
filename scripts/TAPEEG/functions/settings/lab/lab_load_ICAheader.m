function header = lab_load_ICAheader(Filename)

try %#ok<TRYNC>
    load(Filename);
end

if exist('header','var') & exist('W','var')
    header.W = W;
else
    header = [];
end
if exist('ICAchans','var')
    header.ICAchans = ICAchans;
elseif isfield(header,'goodchans')
    header.ICAchans = header.goodchans;
end
if isfield(header,'ICAchans') & isfield(header,'channels')
    header.nactivations = length(header.ICAchans);
end

[~,Filepath,~,FilenameS] = lab_filename(Filename);
if exist(fullfile(Filepath,[FilenameS '_exclude.txt']),'file') & strcmp(FilenameS(end-2:end),'ICA')
    fid = fopen(fullfile(Filepath,[FilenameS '_exclude.txt']));
    badactivations = textscan(fid,'%f');
    badactivations = badactivations{1}(:)';
    fclose(fid);
    clearvars fid
    if ~isempty(badactivations)
        disp(['     read bad activations from ' fullfile(Filepath,[FilenameS '_exclude.txt'])])
        header.badchans = badactivations;
        header.goodchans = setdiff(1:header.numdatachannels,header.badchans);
        header.goodchans = header.goodchans(:)';
    end
elseif exist(fullfile(Filepath,[FilenameS '.info']),'file')
    EEGinfo = lab_read_eeginfo(fullfile(Filepath,[FilenameS '.info']));
    if isfield(EEGinfo,'badchans') & ~isempty(EEGinfo.badchans)
        header.badchans = EEGinfo.badchans;
        header.goodchans = setdiff(1:header.numdatachannels,header.badchans);
        header.goodchans = header.goodchans(:)';
    end
end