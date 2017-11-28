function [cfg,skipprocessing] = lab_set_scale_data(cfg,header)

skipprocessing = 0;

global THeader
if ~isempty(THeader)
    header = THeader;
elseif ~exist('header','var')
    header = [];
end

if ~exist('cfg','var')
    cfg = [];
end

disp ('   Ask for scale settings')
titles={'From','To','Scale factor'};
if exist('cfg','var') & isfield(cfg,'SCALE') & isfield(cfg.SCALE,'scales')
    table = cfg.SCALE.scales;
elseif exist('header','var') & isfield(header,'datatype') & strcmp(header.datatype,'meg') & header.numdatachannels == 306
    if isfield(header,'ecg_ch') & header.ecg_ch > 0
        table = [1,204,1000000000000;205,306,10000000000000;307,header.ecg_ch,10000];
    else
        table = [1,204,1000000000000;205,306,10000000000000];
    end
elseif exist('header','var') & isfield(header,'numdatachannels')
    table = [1,header.numdatachannels,1];
    if isfield(header,'numauxchannels') & header.numauxchannels > 0
        table = cat(1,table,[header.numdatachannels+1,header.numchannels,1]);
    end
elseif exist('cfg','var') & isfield(cfg,'numdatachans')
    if cfg.numdatachans == 306
        table = [1,204,1000000000000;205,306,10000000000000];
    else
        table = [1,cfg.numdatachans,1];
    end
elseif exist('cfg','var') & isfield(cfg,'EXTRA') & isfield(cfg.EXTRA,'numdatachans')
    if cfg.EXTRA.numdatachans == 306
        table = [1,204,1000000000000;205,306,10000000000000];
    else
        table = [1,cfg.EXTRA.numdatachans,1];
    end
else
    table = [1,1,1];
end
tinfo = lab_table_dialog(table,titles,'Scale settings',1);
cfg.SCALE.scales = tinfo;
clearvars info table titles
pause(0.2);
if isempty(cfg.SCALE.scales)
    cfg = rmfield(cfg,'SCALE');
    skipprocessing = 1;
end