function [cfg,skipprocessing] = lab_set_write_result(cfg)

skipprocessing = 0;

if ~exist('cfg','var') | ~isfield(cfg,'RES') | ~isfield(cfg.RES,'format')
    cfg.RES.folder = 'Result';
    cfg.RES.format = 'sef';
    cfg.RES.interpolatemethod = 'spherical';
    cfg.RES.exclude = true;
    cfg.RES.splitchans = 0;
    cfg.RES.appendix = '';
end

Prompt = {};
Formats = {};

Prompt(end+1,:) = {'Output-folder', 'folder'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'text';
Formats(end,1).size = 100;
Formats(end,1).span = [1 2];

Prompt(end+1,:) = {'Interpolate method','interpolatemethod'};
Formats(end+1,1).type = 'list';
Formats(end,1).style = 'popupmenu';
Formats(end,1).format = 'input';
Formats(end,1).items = {'disabled';'spherical';'3D'};
Formats(end,1).span = [1 2];

Prompt(end+1,:) = {'',[]};
Formats(end+1,1).type = 'text';
Formats(end,1).span = [1 2];

Prompt(end+1,:) = {'Write only included channels' 'exclude'};
Formats(end+1,1).type = 'check';
Formats(end,1).size = [-1 -1];
Formats(end,1).span = [1 2];

Prompt(end+1,:) = {'',[]};
Formats(end+1,1).type = 'text';
Formats(end,1).span = [1 2];

Prompt(end+1,:) = {'File format','format'};
Formats(end+1,1).type = 'list';
Formats(end,1).style = 'listbox';
Formats(end,1).format = 'input';
Formats(end,1).items = {'sef';'edf';'eph';'ep';'txt';'fiff'};
Formats(end,1).limits = [0 2]; % multi-select
Formats(end,1).size = [60 100];
Formats(end,1).callback = {@lab_get_format,'scaletxt','format','scaletxt'};

Prompt(end+1,:) = {'Scale TXT','scaletxt'};
Formats(end+1,1).type = 'check';
Formats(end,1).callback = {@lab_show_result,'scaletxt','scaletxt'};

Prompt(end+1,:) = {'',[]};
Formats(end+1,1).type = 'text';
Formats(end,1).span = [1 2];

Prompt(end+1,:) = {'Split output in blocks of ... channels', 'splitchans'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'integer';
Formats(end,1).limits = [0 99];
Formats(end,1).size = 40;
Formats(end,1).span = [1 2];

Prompt(end+1,:) = {'Appendix', 'appendix'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'text';
Formats(end,1).size = 60;
Formats(end,1).span = [1 2];

[cfg.RES,Cancelled] = inputsdlg(Prompt,'Write result',Formats,cfg.RES);
if isempty(cfg.RES) | Cancelled == 1
    cfg.RES = [];
    skipprocessing = 1;
    return
end

return