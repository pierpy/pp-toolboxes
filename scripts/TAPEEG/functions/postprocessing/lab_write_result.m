% Store eeg/meg in different file formats
%
% [cfg] = lab_write_result(data,header,cfg)
%
% data     = matrix (chans x timeframes)
% header   = output of lab_read_data
% cfg      = structure with config (optional)
%
% written by F. Hatz 2012

function [cfg] = lab_write_result(data,header,cfg)

if ~exist('header','var') | ~isfield(header,'samplingrate')
    header = lab_create_header(data);
end

if ~isfield(cfg,'Output_file') & isfield(cfg,'EEG_file')
    cfg.Output_file = cfg.EEG_file;
    cfg.Output_filepath = cfg.EEG_filepath;
else
    cfg.Output_file = header.EEG_file;
    cfg.Output_filepath = header.EEG_filepath;
end
[~,~,~,cfg.Output_fileS] = lab_filename(cfg.Output_file);

if isfield(cfg,'SKIP') & isfield(cfg.SKIP,'RES') & cfg.SKIP.RES == true;
    return
end

disp('Write result')
if ~exist('cfg','var') || ~isfield(cfg,'RES')
    [cfg,skipprocessing] = lab_set_write_result(cfg);
    if skipprocessing == 1
        return
    end
end

if ~isfield(cfg.RES,'appendix')
    cfg.RES.appendix = '_result';
end

% Interpolate
if ~strcmp(cfg.RES.interpolatemethod,'disabled')
    [data,header] = lab_interpolate_bad(data,header,cfg.RES.interpolatemethod);
end

% Eliminate cfg.exclude
if cfg.RES.exclude == 1 & isfield(cfg,'exclude') & cfg.exclude > 0
    [data,header] = lab_reduce_channels(data,header,setdiff(1:header.numchannels,cfg.exclude));
end

Output_filepath = cfg.Output_filepath;
if isfield(cfg.RES,'folder') & ~isempty(cfg.RES.folder)
    Output_filepath = fullfile(Output_filepath,cfg.RES.folder);
    warning off %#ok<WNOFF>
    mkdir(Output_filepath);
    warning on %#ok<WNON>
end

% Separate in blocks
if cfg.RES.splitchans > 0
    if cfg.RES.splitchans == 102 & strcmp(header.datatype,'meg') & header.numdatachannels == 306
        blocks = [1 102;103 204;205 306];
        blocknames = {'transversal';'radial';'mag'};
    else
        tmp = ceil(header.numdatachannels / cfg.RES.splitchans);
        for i = 1:tmp
            blocks(i,1) = (i-1)*cfg.RES.splitchans + 1;
            blocks(i,2) = i*cfg.RES.splitchans;
            blocknames{i,1} = ['block' num2str(i)];
        end
        blocks(end,2) = header.numdatachannels;
        clearvars tmp
    end
    for i = 1:size(blocks,1)
        if header.numchannels > header.numdatachannels
            chans = union((blocks(i,1):blocks(i,2)),(header.numdatachannels+1:header.numchannels));
        else
            chans = (blocks(i,1):blocks(i,2));
        end
        Output_file=[cfg.Output_fileS '_' blocknames{i,1} cfg.RES.appendix '.sef'];
        [data2,header2] = lab_reduce_channels(data,header,chans);
        for j = 1:length(cfg.RES.format)
            [filename,cfg.RES] = lab_save_data(data2,header2,cfg.RES.format{j},fullfile(Output_filepath,Output_file),cfg.RES,cfg);
        end
    end
else
    Output_file=[cfg.Output_fileS cfg.RES.appendix '.sef'];
    for i = 1:length(cfg.RES.format)
        [filename,cfg.RES] = lab_save_data(data,header,cfg.RES.format{i},fullfile(Output_filepath,Output_file),cfg.RES,cfg);
    end
end

%--------------------------------------------------------------------------
% Write verbose file (*.vrb)
%--------------------------------------------------------------------------
fid=fopen(fullfile(Output_filepath,[cfg.Output_fileS '_writeresult.vrb']),'w');
fprintf(fid,'Write result\n');
fprintf(fid,['EEG-file: ' cfg.Output_file]);
fprintf(fid,'\n\n');
fprintf(fid,['File format: ' sprintf('%s|',cfg.RES.format{:})]);
fprintf(fid,'\n');
fprintf(fid,['Interpolate bad (method): ' cfg.RES.interpolatemethod]);
fprintf(fid,'\n');
if cfg.RES.exclude == true
    fprintf(fid,'Write only included channels: yes');
else
    fprintf(fid,'Write only included channels: no');
end
fprintf(fid,'\n');
fprintf(fid,['Appendix for files: ' cfg.RES.appendix]);
fprintf(fid,'\n');
if ~isempty(cfg.RES.scaletxt)
    fprintf(fid,['Scale factor for txt: ' cfg.RES.scaletxt]);
    fprintf(fid,'\n');
end
fclose(fid);

cfg.SKIP.RES = true;

return