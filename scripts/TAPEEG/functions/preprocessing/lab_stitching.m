% Script for stitching of eeg/meg-files
%
% [data,header,cfg,skipprocessing] = lab_stitching(data,header,cfg,calc)
%
% written by F. Hatz 2013

function [data,header,cfg,skipprocessing] = lab_stitching(data,header,cfg,calc)

skipprocessing = 0;
global stitchdata stitchheader stitchinfo

if ~exist('cfg','var')
    cfg = [];
end

if ~isfield(cfg,'STITCH') | ~isfield(cfg.STITCH,'length')
    [cfg,skipprocessing] = lab_set_stitching(cfg);
    if skipprocessing == 1
        return
    end
end

if cfg.STITCH.length == 0 & cfg.STITCH.stitchall == false
    % exist without stitching
    return
end

if isempty(stitchinfo) | isempty(find(strcmp(stitchinfo(:,1),fullfile(cfg.EEG_filepath,cfg.EEG_file)),1,'first'))
    lab_prepare_stitching(cfg,calc);
end
findex = find(strcmp(stitchinfo(:,1),fullfile(cfg.EEG_filepath,cfg.EEG_file)),1,'first');
if ~isempty(findex) & stitchinfo{findex,2} == 1
    dostitch = 1;
else
    dostitch = 0;
end

% stiching of short segments (if enabled)
if ~isempty(stitchdata) & ~isempty(stitchheader)
    [data,header] = lab_stitch(stitchdata,stitchheader,data,header,cfg.STITCH.window);
    stitchdata = [];
    stitchheader = [];
end
if dostitch == 1
    stitchdata = data;
    stitchheader = header;
    skipprocessing = 1;
end

% generate new name
if isfield(cfg,'patient') & ~isempty(cfg.patient)
    [~,~,format,fileS] = lab_filename(header.EEG_file);
    if cfg.STITCH.stitchall == false
        header.EEG_file = [fileS '_Stitch' stitchinfo{findex,3} '.' format];
    else
        header.EEG_file = [fileS '_Stitch.' format];
    end
    cfg.EEG_file = header.EEG_file;
    cfg.Output_file = header.EEG_file;
end