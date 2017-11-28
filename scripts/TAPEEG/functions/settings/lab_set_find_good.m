function [cfg,skipprocessing] = lab_set_find_good(cfg,header)

global THeader
if ~isempty(THeader)
    header = THeader;
elseif ~exist('header','var')
    header = [];
end
    
if ~exist('cfg','var')
    cfg = [];
end

if ~isfield(cfg,'BADELEC')
    cfg.BADELEC = [];
end
if isempty(cfg.BADELEC)
    cfg.BADELEC.filemethod = 'skip';
    cfg.BADELEC.length = 4; % set length of epochs to 4 seconds
    cfg.BADELEC.percentbad = 70;
    cfg.BADELEC.freqlim50 = 50;
    cfg.BADELEC.freqlim60 = 50;
    cfg.BADELEC.freqlimlow = 70;
    cfg.BADELEC.freqlimhigh = 50;
    cfg.BADELEC.spectslow = [0.5 2];
    cfg.BADELEC.spectshigh = [15 50];
    cfg.BADELEC.zvaluebroken = 4;
    if isfield(header,'interpolated') & ~isempty(header.interpolated)
        cfg.BADELEC.zvaluevars = 4;
    else
        cfg.BADELEC.zvaluevars = 3;
    end
    cfg.BADELEC.zvaluehurst = 3;
    cfg.BADELEC.zvaluemedian = 0;
    cfg.BADELEC.zvaluekurtosis = 0;
    if ~isfield(header,'numdatachannels') | header.numdatachannels > 72
        cfg.BADELEC.zvaluecorr = 3;
    else
        cfg.BADELEC.zvaluecorr = 0;
    end
    cfg.BADELEC.LAPL.lap_maxdistance = 2.5;
    cfg.BADELEC.LAPL.lap_weightmaxdistance = 30;
    cfg.BADELEC.LAPL.lap_excluderef = true;
end
[cfg.BADELEC,skipprocessing] = lab_set_detect_bad(cfg.BADELEC,header,cfg,0,0,2,0,1,0,1,1,0,1);