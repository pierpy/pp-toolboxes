% Coregistration of MRI-file and electrodes/sensors information
%
% [LOCS,cfg] = lab_coregistration(data,header,cfg)
%
% Written by F. Hatz 2013

function [LOCS,cfg] = lab_coregistration(~,header,cfg)

if ~exist('cfg','var')
    cfg = [];
end
if ~exist('header','var') | ~isfield(header,'locs')
    disp('Abort: no location information')
    return
end
LOCS = header.locs;

if isfield(cfg,'SKIP') & isfield(cfg.SKIP,'COREG') & cfg.SKIP.COREG == true;
    return
end
if ~isfield(cfg,'COREG') | isempty(cfg.COREG) | ~isfield(cfg.COREG,'MRI_file')
    [cfg.COREG,skipprocessing] = lab_set_coreg([],header,1,0);
    if skipprocessing == 1
        LOCS = [];
        return
    else
        pause(0.2);
    end
end

if ~isfield(header,'EEG_filepath') & isfield(cfg,'Output_filepath')
    header.EEG_filepath = cfg.Output_filepath;
    header.EEG_file = cfg.Output_file;
end
if ~isfield(header,'EEG_filepath') | ~exist(header.EEG_filepath,'dir')
    disp('Abort: no FilePath provided')
    return
end

tmp = strfind(cfg.COREG.MRI_file,filesep);
if isempty(tmp) & ~isempty(cfg.COREG.MRI_file)
    findmri = dir(fullfile(header.EEG_filepath,['*.' cfg.COREG.MRI_file]));
    if size(findmri,1) > 0
        MRI_file = fullfile(header.EEG_filepath,findmri(1,1).name);
    else
        MRI_file = '';
    end
elseif exist(cfg.COREG.MRI_file,'file')
    MRI_file = cfg.COREG.MRI_file;
else
    MRI_file = '';
end
if isempty(MRI_file) | ~exist(MRI_file,'file')
    disp('Abort: no valid MRI-file found')
    return
end

[LOCS,~,cfg.COREG] = lab_coreg(MRI_file,header.locs,cfg.COREG);