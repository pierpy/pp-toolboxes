function [cfg,skipprocessing] = lab_set_stitchingall(cfg)

skipprocessing = 0;

if ~exist('cfg','var') | ~isfield(cfg,'STITCHALL') | ~isfield(cfg.STITCHALL,'window')
    cfg.STITCHALL.window = 1;
    cfg.STITCHALL.enable = true;
    cfg.STITCHALL.interpolatebad = 'spherical';
end

Prompt = cell(0,2);
Formats = {};

Prompt(end+1,:) = {'Stitch all files in folder' 'enable'};
Formats(end+1,1).type = 'check';

Prompt(end+1,:) = {'overlap window (seconds)','window'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'float';
Formats(end,1).limits = [0 inf];
Formats(end,1).size = 30;

Prompt(end+1,:) = {'Interpolate bad channels' 'interpolatebad'};
Formats(end+1,1).type = 'list';
Formats(end,1).style = 'popupmenu';
Formats(end,1).format = 'input';
Formats(end,1).items = {'disabled';'spherical';'3D'};
Formats(end,1).span = [1 2];

[cfg.STITCHALL,Cancelled] = inputsdlg(Prompt,'Stitch folder setting',Formats,cfg.STITCHALL,2);
if Cancelled == 1
    cfg.STITCHALL = [];
    skipprocessing = 1;
    return
end
if cfg.STITCHALL.enable == false
    cfg.STITCHALL = [];
end

return