% Script for stitching of eeg/meg-files
%
% [data,header,cfg,skipprocessing] = lab_stitchingall(data,header,cfg)
%
% written by F. Hatz 2013

function [data,header,cfg,skipprocessing] = lab_stitchingall(data,header,cfg)

global stitchdataALL stitchheaderALL

if ~exist('skipprocessing','var')
    skipprocessing = 0;
end

if ~exist('cfg','var')
    cfg = [];
end

if ~isfield(cfg,'STITCHALL') | ~isfield(cfg.STITCHALL,'enable')
    [cfg,skipprocessing] = lab_set_stitchingall(cfg);
    if skipprocessing == 1
        return
    end
end

if cfg.STITCHALL.enable == false
    % exist without stitching
    return
else
    disp('Stitching All')
end

if isfield(cfg.STITCHALL,'interpolatebad') & ~strcmp(cfg.STITCHALL.interpolatebad,'disabled') & ~isempty(data)
    [data,header] = lab_interpolate_bad(data,header,cfg.STITCHALL.interpolatebad);
end

% stiching of segments
if ~isempty(stitchdataALL) & ~isempty(stitchheaderALL) & ~isempty(data)
    [stitchdataALL,stitchheaderALL] = lab_stitch(stitchdataALL,stitchheaderALL,data,header,cfg.STITCHALL.window);
elseif ~isempty(data)
    stitchdataALL = data;
    stitchheaderALL = header;
end

if ~isfield(cfg,'lastsegment') | cfg.lastsegment == true
    data = stitchdataALL;
    header = stitchheaderALL;
    stitchdataALL = [];
    stitchheaderALL = [];
else
    disp('    Keep data in memory and proceed with next segment/file')
    data = [];
    header = [];
    skipprocessing = 1;
    return
end

% generate new name
if isfield(cfg,'patient') & ~isempty(cfg.patient)
    [~,~,format,fileS] = lab_filename(cfg.EEG_file);
    tmp = strfind(fileS,'_Stitch');
    if ~isempty(tmp)
        fileS = fileS(1:tmp(1)-1);
    end
    clearvars tmp
    cfg.EEG_file = [fileS '_Stitch.' format];
    header.EEG_file = cfg.EEG_file;
    cfg.Output_file = header.EEG_file;
end