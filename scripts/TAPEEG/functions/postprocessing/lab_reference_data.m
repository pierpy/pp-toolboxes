% Calculate new references and store results
%
% [dataR,headerR,cfg] = lab_reference_data(data,header,cfg)
%
% data    = matrix (chans x timeframes)
% header  = output of lab_read_data
% cfg     = structure with config (optional)
%
% written by F. Hatz 2012

function [dataR,headerR,cfg] = lab_reference_data(data,header,cfg)

dataR = [];
headerR = [];
if isempty(data)
    return
end
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

if isfield(cfg,'SKIP') & isfield(cfg.SKIP,'REF') & cfg.SKIP.REF == true;
    return
end

disp('Re-referencing EEG-data')
if ~isfield(cfg,'REF') ||  ~isfield(cfg.REF,'format')
    [cfg,skipprocessing]  = lab_set_reference_data(cfg,header);
    if skipprocessing == 1
        return
    end
end

%--------------------------------------------------------------------------
% Correct for missing good/bad channels
%--------------------------------------------------------------------------
if ~isfield(header,'goodchans')
    header.goodchans = 1:header.numdatachannels;
    header.badchans = [];
end

%--------------------------------------------------------------------------
% Correct bad/good ref_chan
%--------------------------------------------------------------------------
if isnumeric(header.ref_chan) & ~isempty(header.ref_chan) & header.ref_chan > 0
    header.goodchans = union(header.goodchans,header.ref_chan);
    header.goodchans = header.goodchans(:)';
    header.badchans = setdiff(header.badchans,header.ref_chan);
    header.badchans = header.badchans(:)';
end

Output_file = cfg.Output_file;
Output_fileS = cfg.Output_fileS;
Output_filepath = cfg.Output_filepath;
if isfield(cfg.REF,'folder') & ~isempty(cfg.REF.folder)
    cfg.Output_filepath = fullfile(cfg.Output_filepath,cfg.REF.folder);
    warning off %#ok<WNOFF>
    mkdir(cfg.Output_filepath);
    warning on %#ok<WNON>
end

cfg.SKIP.REF = true;
SKIP = cfg.SKIP;
SKIPtmp = cfg.SKIP;

%--------------------------------------------------------------------------
% Write corrected EEG mean
%--------------------------------------------------------------------------
if max(strcmp(cfg.REF.eegsource,'mean'))
    [dataR,headerR] = lab_references(data,header,'mean');
    if cfg.REF.interpolatebad == 1
        [dataR,headerR] = lab_interpolate_bad(dataR,headerR);
    end
    disp('   Write mean-file')
    filename = fullfile(cfg.Output_filepath,[Output_fileS '_mean.sef']);
    for i = 1:length(cfg.REF.format)
        [cfg.Output_file,cfg.REF] = lab_save_data(dataR,headerR,cfg.REF.format{i},filename,cfg.REF,cfg);
        cfg.Output_file = lab_filename(cfg.Output_file);
    end
    cfg.SKIP = SKIPtmp;
    cfg = lab_postprocessing(dataR,headerR,cfg,true);
    SKIP = merge_skip(SKIP,cfg.SKIP);
end
cfg.Output_file=Output_file;
cfg.Output_fileS=Output_fileS;

%--------------------------------------------------------------------------
% Write corrected EEG median
%--------------------------------------------------------------------------
if max(strcmp(cfg.REF.eegsource,'median'))
    [dataR,headerR] = lab_references(data,header,'median');
    if cfg.REF.interpolatebad == 1
        [dataR,headerR] = lab_interpolate_bad(dataR,headerR);
    end
    disp('   Write median-file')
    filename = fullfile(cfg.Output_filepath,[Output_fileS '_median.sef']);
    for i = 1:length(cfg.REF.format)
        [cfg.Output_file,cfg.REF] = lab_save_data(dataR,headerR,cfg.REF.format{i},filename,cfg.REF,cfg);
        cfg.Output_file = lab_filename(cfg.Output_file);
    end
    cfg.SKIP = SKIPtmp;
    cfg = lab_postprocessing(dataR,headerR,cfg,true);
    SKIP = merge_skip(SKIP,cfg.SKIP);
end
cfg.Output_file=Output_file;
cfg.Output_fileS=Output_fileS;

%--------------------------------------------------------------------------
% Write corrected EEG laplacian
%--------------------------------------------------------------------------
if max(strcmp(cfg.REF.eegsource,'laplacian'))
    [dataR,headerR,cfg.REF] = lab_references(data,header,'laplacian',cfg.REF);
    if cfg.REF.interpolatebad == 1
        [dataR,headerR] = lab_interpolate_bad(dataR,headerR);
    end
    disp('   Write laplacian-file')
    filename = fullfile(cfg.Output_filepath,[Output_fileS '_laplacian.sef']);
    for i = 1:length(cfg.REF.format)
        [cfg.Output_file,cfg.REF] = lab_save_data(dataR,headerR,cfg.REF.format{i},filename,cfg.REF,cfg);
        cfg.Output_file = lab_filename(cfg.Output_file);
    end
    cfg.SKIP = SKIPtmp;
    cfg = lab_postprocessing(dataR,headerR,cfg,true);
    SKIP = merge_skip(SKIP,cfg.SKIP);
end
cfg.Output_file=Output_file;
cfg.Output_fileS=Output_fileS;

%--------------------------------------------------------------------------
% Write corrected EEG montage
%--------------------------------------------------------------------------
if max(strcmp(cfg.REF.eegsource,'montage'))
    % Correct montage
    if ~isfield(cfg.REF,'montage') | isempty(cfg.REF.montage)
        disp('   no montage loaded, use input structure')
        cfg.REF.montage = lab_create_montage(size(data,1),header);
    end
    if cfg.REF.montage(1,1).numchans ~= header.numchannels
        cfg.REF.montage = lab_reduce_montage(cfg.REF.montage,cfg,header,true);
    end
    if cfg.REF.montage(1,1).numchans > header.numchannels
        montage = lab_create_montage(size(data,1),header);
    else
        montage = cfg.REF.montage;
    end
    if cfg.REF.interpolatebad == 1
        [dataR,headerR] = lab_interpolate_bad(data,header);
    else
        dataR = data;
        headerR = header;
    end
    disp('   Write montage-file')
    for i = 1:size(montage,2);
        [dataM,headerM] = lab_references(dataR,headerR,montage(1,i));
        filename = fullfile(cfg.Output_filepath,[Output_fileS '_' montage(1,i).name '.sef']);
        for j = 1:length(cfg.REF.format)
            [cfg.Output_file,cfg.REF] = lab_save_data(dataM,headerM,cfg.REF.format{j},filename,cfg.REF,cfg);
            cfg.Output_file = lab_filename(cfg.Output_file);
        end
        cfg.SKIP = SKIPtmp;
        cfg = lab_postprocessing(dataR,headerR,cfg,true);
        SKIP = merge_skip(SKIP,cfg.SKIP);
    end
    dataR = dataM;
    headerR = headerM;
    clearvars dataM headerM
end
cfg.Output_file=Output_file;
cfg.Output_fileS=Output_fileS;

%--------------------------------------------------------------------------
% Write corrected EEG channel reference
%--------------------------------------------------------------------------
for j = 1:length(cfg.REF.eegsource)
    if isnumeric(cfg.REF.eegsource{j}) | ~isnan(str2double(cfg.REF.eegsource{j}))
        [dataR,headerR,cfg.REF] = lab_references(data,header,cfg.REF.eegsource{j},cfg.REF);
        if cfg.REF.interpolatebad == 1
            [dataR,headerR] = lab_interpolate_bad(dataR,headerR);
        end
        disp('   Write channel reference-file')
        filename = fullfile(cfg.Output_filepath,[Output_fileS '_chanref.sef']);
        for i = 1:length(cfg.REF.format)
            [cfg.Output_file,cfg.REF] = lab_save_data(dataR,headerR,cfg.REF.format{i},filename,cfg.REF,cfg);
            cfg.Output_file = lab_filename(cfg.Output_file);
        end
        cfg.SKIP = SKIPtmp;
        cfg = lab_postprocessing(dataR,headerR,cfg,true);
        SKIP = merge_skip(SKIP,cfg.SKIP);
    end
end
cfg.Output_file=Output_file;
cfg.Output_fileS=Output_fileS;
cfg.Output_filepath = Output_filepath;
cfg.SKIP = SKIP;

%--------------------------------------------------------------------------
% Write verbose file (*.vrb)
%--------------------------------------------------------------------------
for i = 1:length(cfg.REF.eegsource)
    if isnumeric(cfg.REF.eegsource{i})
        cfg.REF.eegsource{i} = num2str(cfg.REF.eegsource{i});
    end
end
fid=fopen(fullfile(cfg.Output_filepath,[cfg.Output_fileS '_newref.vrb']),'w');
fprintf(fid,'Save new references\n');
fprintf(fid,['EEG-file: ' cfg.Output_file]);
fprintf(fid,'\n\n');
fprintf(fid,['New references: ' sprintf('%s|',cfg.REF.eegsource{:})]);
fprintf(fid,'\n\n');
if cfg.REF.interpolatebad == 1
    fprintf(fid,'Interpolate bad channels: enabled');
else
    fprintf(fid,'Interpolate bad channels: disabled');
end
fprintf(fid,'\n\n');
if isfield(cfg.REF,'LAPL') & isfield(cfg.REF.LAPL,'lap_maxdistance')
    fprintf(fid,['Laplacian Maxdistance (number of minimal electrode distances): ' num2str(cfg.REF.LAPL.lap_maxdistance)]);
    fprintf(fid,'\n');
    fprintf(fid,['Laplacian Weight Maxdistance (percent): ' num2str(cfg.REF.LAPL.lap_weightmaxdistance)]);
    fprintf(fid,'\n');
end
fclose(fid);

end

function SKIP = merge_skip(SKIP,SKIP2)

if ~isstruct(SKIP) | ~isstruct(SKIP2) | isempty(SKIP2)
    return
end
    
Fields = fieldnames(SKIP2);
for i = 1:length(Fields)
    if SKIP2.(Fields{i}) == true
        SKIP.(Fields{i}) = true;
    end
end

end