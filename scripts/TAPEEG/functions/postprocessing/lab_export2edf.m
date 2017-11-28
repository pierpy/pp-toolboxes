% Store EEG/MEG as EDF
%
% cfg = lab_export2edf(data,header,cfg,nofolder)
%
% data     = eeg/meg data (chans x timeframes)
% header   = output of lab_read_data
% cfg      = structure with config (optional)
% nofolder = if variable is defined, config for specific folder to store
%            EDF is ignored
%
% Written by F. Hatz 2012 Neurology Basel

function [data_edf,header_edf,cfg] = lab_export2edf(data,header,cfg,nofolder)

if ~exist('header','var') | ~isfield(header,'samplingrate')
    header = lab_create_header(data);
end
if ~exist('cfg','var') || ~isfield(cfg,'settings_path')
    cfg.settings_path = pwd;
end

if ~isfield(cfg,'Output_file')
    cfg.Output_file = header.EEG_file;
    cfg.Output_filepath = header.EEG_filepath;
end
[~,~,~,cfg.Output_fileS] = lab_filename(cfg.Output_file);

if ~isfield(cfg,'EDF_file')
    cfg.EDF_file = [cfg.Output_fileS '.edf'];
    cfg.EDF_filepath = cfg.Output_filepath;
end

if isfield(cfg,'SKIP') & isfield(cfg.SKIP,'EDF') & cfg.SKIP.EDF == true;
    data_edf = [];
    header_edf = [];
    return
end

disp('Export to EDF-file')
if ~isfield(cfg,'EDF')
    [cfg,skipprocessing] = lab_set_export2edf(cfg,header);
    if skipprocessing == 1
        data_edf = [];
        header_edf = [];
        return
    end
end

% Calculate Reference
if ~isfield(header,'numdatachannels')
    header.numdatachannels = header.numchannels;
end
if strcmp(cfg.EDF.eegsource,'montage') & isfield(cfg.EDF,'montage')
    % Interpolate bad channels
    if cfg.EDF.interpolatebad == 1
        data = lab_interpolate_bad(data,header);
    end
    % Correct montage
    if ~isfield(cfg.EDF,'montage') | isempty(cfg.EDF.montage)
        disp('   no montage loaded, use input structure')
        cfg.EDF.montage = lab_create_montage(size(data,1),header);
    end
    if cfg.EDF.montage(1,1).numchans ~= header.numchannels
        cfg.EDF.montage = lab_reduce_montage(cfg.EDF.montage,cfg,header,true);
    end
    if cfg.EDF.montage(1,1).numchans > header.numchannels
        montage = lab_create_montage(size(data,1),header);
    else
        montage = cfg.EDF.montage;
    end
    % Calculate montage
    if montage(1,1).numchans <= header.numchannels
        [data_edf,header_edf] = lab_references(data,header,montage(1,1));
    else
        disp('    Montage file not matching, using signal reference')
        data_edf = data;
        header_edf = header;
        cfg.EDF.eegsource = 'input';
    end
else
    % reference data
    [data_edf,header_edf] = lab_references(data,header,cfg.EDF.eegsource);
    % Interpolate bad channels
    if cfg.EDF.interpolatebad == 1
        data_edf = lab_interpolate_bad(data_edf,header_edf);
    end
    % Reduce channels to 255
    if size(data,1) > 255
        disp('    For edf-export channels are reduced to first 255 channels')
        [data_edf,header_edf] = lab_reduce_channels(data_edf,header_edf,1:255);
    end
end

% Write EDF-File
disp('   Write EDF-file')
header_edf.numchannels = size(data_edf,1);
header_edf.numtimeframes = size(data_edf,2);
if ~isfield(header_edf,'subject')
    header_edf.subject.name = cfg.Output_fileS;
    header_edf.subject.name = regexprep(header_edf.subject.name,' ','_');
end
if size(header_edf.channels,1) > size(data_edf,1)
    if strcmp(header_edf.channels(end,1:3),'ECG') && isfield(header,'ecg_ch') && header.ecg_ch > 0
        data_edf(end+1,:) = data(header.ecg_ch,:);
        header_edf.numdatachannels = header_edf.numchannels;
        header_edf.numauxchannels = 1;
        header_edf.numchannels = size(data_edf,1);
        header_edf.ecg_ch = size(data_edf,1);
    end
    header_edf.channels = header_edf.channels(1:size(data_edf,1),1);
end

if isfield(cfg.EDF,'mountpath') & ~isempty(cfg.EDF.mountpath)
    system(cfg.EDF.mountpath);
end
if ~isempty(cfg.EDF.filepath) & exist(cfg.EDF.filepath,'dir')
    cd(cfg.EDF.filepath);
    warning off %#ok<WNOFF>
    mkdir(cfg.EDF_file(1:end-4));
    warning on %#ok<WNON>
    tmp = strfind(cfg.EDF.filepath,filesep);
    if tmp(end) == length(cfg.EDF.filepath);
        cfg.EDF_filepath = [cfg.EDF.filepath cfg.EDF_file(1:end-4) filesep ];
    else
        cfg.EDF_filepath = [cfg.EDF.filepath filesep cfg.EDF_file(1:end-4) filesep];
    end
    clearvars tmp oldFolder
else
    if ~exist('nofolder','var')
        tmp = strfind(cfg.EDF_filepath,filesep);
        if tmp(end) ~= length(cfg.EDF_filepath);
            cfg.EDF_filepath = [cfg.EDF_filepath filesep];
        end
        warning off %#ok<WNOFF>
        mkdir ([cfg.EDF_filepath 'EDF']);
        warning on %#ok<WNON>
        cfg.EDF_filepath = [cfg.EDF_filepath 'EDF' filesep];
    end
end
lab_write_edf(fullfile(cfg.EDF_filepath,cfg.EDF_file),data_edf,header_edf);

if ~strcmp(cfg.EDF_file(end-4),'~')
    %--------------------------------------------------------------------------
    % Write verbose file (*.vrb)
    %--------------------------------------------------------------------------
    fid=fopen(fullfile(cfg.EDF_filepath,[cfg.EDF_file(1:end-4) '_EDF.vrb']),'w');
    fprintf(fid,'EDF-file\n');
    fprintf(fid,'\n');
    if isnumeric(cfg.EDF.eegsource)
        fprintf(fid,['EEG reference: ' num2str(cfg.EDF.eegsource)]);
    else
        fprintf(fid,['EEG reference: ' cfg.EDF.eegsource]);
    end
    fprintf(fid,'\n\n');
    if strcmp(cfg.EDF.eegsource,'montage') & isfield(cfg.EDF,'montage')
        fprintf(fid,lab_montage2txt(cfg.EDF.montage(1,1),header));
    end
end

cfg.SKIP.EDF = true;

end