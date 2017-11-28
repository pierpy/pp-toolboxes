function [cfg,skipprocessing] = lab_set_average(cfg,header,data,short)

skipprocessing = 0;

if ~exist('short','var')
    short = 0;
end
if ~exist('data','var')
    data = [];
end

global THeader TData
if ~isempty(TData)
    data = TData;
elseif ~exist('data','var')
    data = [];
end
if ~isempty(THeader)
    header = THeader;
elseif ~exist('header','var')
    header = [];
end

if ~exist('cfg','var') | isempty(cfg) | ~isfield(cfg,'AVG') | ~isfield(cfg.AVG,'marker')
    cfg.AVG.marker = {'*'};
    cfg.AVG.combinemarker = true;
    cfg.AVG.markerlength = 1000;
    cfg.AVG.markerOffset = 0;
    cfg.AVG.correctbaseline = 'disabled';
    cfg.AVG.AVGmethod = 'mean';
    cfg.AVG.Reject.routine = 'FASTER';
    cfg.AVG.Reject.BAD = 1.5;
    if short == 0
        cfg.AVG.Reject.method = 'percent';
        cfg.AVG.Reject.percent = 5;
    else
        cfg.AVG.Reject.method = 'single';
        cfg.AVG.Reject.percent = [];
    end
    if short == 0
        cfg.AVG.folder = 'AVG';
        cfg.AVG.format = {'eph'};
        cfg.AVG.scaletxt = [];
        cfg.AVG.interpolate = true;
        cfg.AVG.Correct = [];
        cfg.AVG.eegsource = {'input'};
        cfg.AVG.FILT = [];
        cfg.AVG.JITTER = [];
    end
end

if isfield(header,'events') & isfield(header.events,'TYP') & ~isempty(header.events.TYP)
    tmp = unique(header.events.TYP);
    markerlist = [cellstr('**edit**') tmp(:)'];
else
    markerlist = [];
end
if isfield(cfg,'MARK') & isfield(cfg.MARK,'edit') & ~isempty(cfg.MARK.edit)
    markerlist = union(markerlist,cfg.MARK.edit(:,7));
end
if (isfield(cfg,'STITCH') & ~isempty(cfg.STITCH)) | (isfield(cfg,'STITCHALL') & ~isempty(cfg.STITCHALL))
    markerlist = union(markerlist,cellstr('Stitch'));
end

Prompt = cell(0,2);
Formats = [];

if short == 0
    refitems = {'input','mean','median','laplacian','montage','channels'};
    if ~max(strcmp({'input';'mean';'median';'laplacian';'montage'},cfg.AVG.eegsource))
        refitems = cat(2,refitems,cfg.AVG.eegsource);
    end

    Prompt(end+1,:) = {'Output-folder', 'folder'};
    Formats(end+1,1).type = 'edit';
    Formats(end,1).format = 'text';
    Formats(end,1).size = 100;
    Formats(end,1).span = [1 3];
    
    Prompt(end+1,:) = {'Reference','eegsource'};
    Formats(end+1,1).type = 'list';
    Formats(end,1).style = 'popupmenu';
    Formats(end,1).format = 'input';
    Formats(end,1).items = refitems;
    Formats(end,1).callback = {@lab_get_eegsource,'@ALL','@ALL',cfg,header};
    
    Prompt(end+1,:) = {'Montage','montage'};
    Formats(end+1,1).type = 'check';
    Formats(end,1).callback = {@lab_load_montage,'montage','montage',cfg,header};
    
    Prompt(end+1,:) = {'Laplacian','LAPL'};
    Formats(end+1,1).type = 'check';
    Formats(end,1).callback = {@lab_get_laplacian,'LAPL','LAPL'};
    
    Prompt(end+1,:) = {'Filter','FILT'};
    Formats(end+1,1).type = 'check';
    Formats(end,1).callback = {@set_filtering,'FILT','FILT'};
    Formats(end,1).span = [1 3];
    
    Prompt(end+1,:) = {'Interpolate bad channels','interpolate'};
    Formats(end+1,1).type = 'check';
    Formats(end,1).span = [1 3];
    
    Prompt(end+1,:) = {'',''};
    Formats(end+1,1).type = 'text';
    Formats(end,1).span = [1 3];
end

Prompt(end+1,:) = {'Markers','marker'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'vector';
Formats(end,1).size = 170;
if exist('markerlist','var')
    Formats(end,1).items = markerlist;
end
Formats(end,1).span = [1 2];

Prompt(end+1,:) = {'Combine markers','combinemarker'};
Formats(end+1,1).type = 'check';

Prompt(end+1,:) = {'Length of pre-marker period','markerOffset'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'integer';
Formats(end,1).limits = [0 9999];
Formats(end,1).size = 35;
Formats(end,1).callback = {@set_markerOffset,'correctbaseline','correctbaseline','markerOffset'};
Formats(end,1).span = [1 2];

Prompt(end+1,:) = {'Length of average','markerlength'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'integer';
Formats(end,1).limits = [0 9999];
Formats(end,1).size = 35;

Prompt(end+1,:) = {'',''};
Formats(end+1,1).type = 'text';
Formats(end,1).span = [1 3];

Prompt(end+1,:) = {'Average method','AVGmethod'};
Formats(end+1,1).type = 'list';
Formats(end,1).style = 'popupmenu';
Formats(end,1).format = 'input';
Formats(end,1).items = {'mean','median'};

Prompt(end+1,:) = {'Reject','Reject'};
Formats(end+1,1).type = 'check';
Formats(end,1).callback = {@set_reject,'Reject','Reject',header,cfg};
Formats(end,1).span = [1 2];

Prompt(end+1,:) = {'Baseline-Correction','correctbaseline'};
Formats(end+1,1).type = 'list';
Formats(end,1).style = 'popupmenu';
Formats(end,1).format = 'input';
Formats(end,1).items = {'disabled','all TFs','pre-marker period'};
Formats(end,1).callback = {@set_correctbaseline,'correctbaseline','correctbaseline','markerOffset'};
Formats(end,1).span = [1 3];

if short == 0
    Prompt(end+1,:) = {'',''};
    Formats(end+1,1).type = 'text';
    Formats(end,1).span = [1 3];
    
    Prompt(end+1,:) = {'Calculate Jitter','JITTER'};
    Formats(end+1,1).type = 'check';
    Formats(end,1).callback = {@lab_get_JITTER,'JITTER','JITTER','markerOffset'};
    Formats(end,1).span = [1 3];
    
    Prompt(end+1,:) = {'Detect and interpolate bad channels in average results','Correct'};
    Formats(end+1,1).type = 'check';
    Formats(end,1).span = [1 3];
    Formats(end,1).callback = {@set_correct,'Correct','Correct',header,cfg};
    
    Prompt(end+1,:) = {'',''};
    Formats(end+1,1).type = 'text';
    Formats(end,1).span = [1 3];
    
    Prompt(end+1,:) = {'File format','format'};
    Formats(end+1,1).type = 'list';
    Formats(end,1).style = 'listbox';
    Formats(end,1).format = 'input';
    Formats(end,1).items = {'sef';'edf';'eph';'ep';'txt'};
    Formats(end,1).limits = [0 2]; % multi-select
    Formats(end,1).size = [60 85];
    Formats(end,1).callback = {@lab_get_format,'scaletxt','format','scaletxt'};
    
    Prompt(end+1,:) = {'Scale TXT','scaletxt'};
    Formats(end+1,1).type = 'check';
    Formats(end,1).callback = {@lab_show_result,'scaletxt','scaletxt'};
    Formats(end,1).span = [1 2];
    
    Prompt(end+1,:) = {'Dipol fitting','IS'};
    Formats(end+1,1).type = 'check';
    Formats(end,1).callback = {@lab_set_dipolfit,'@ALL','@ALL',header,data,1};
end

[cfg.AVG,Cancelled] = inputsdlg(Prompt,'Average',Formats,cfg.AVG,2);
if isempty(cfg.AVG) | Cancelled == 1
    cfg.AVG = [];
    skipprocessing = 1;
    return
end

end

function settings = set_filtering(settings)
if isempty(settings)
    settings.filtermode = 'firls';
    settings.filtorder = 0;
    settings.highpass = 1;
    settings.lowpass = 30;
end
[settings,Cancelled] = lab_set_filter(settings,1);
if Cancelled == 1
    settings = [];
end
end

function settings = set_reject(settings,header,cfg)
    if ~isfield(settings,'method')
        settings.routine = 'FASTER';
        settings.BAD = 1.5;
        settings.method = 'percent';
        settings.percent = 5;
    end
    Prompt = cell(0,3);
    Formats = [];
    
    Prompt(end+1,:) = {'Routine','routine',''};
    Formats(end+1,1).type = 'list';
    Formats(end,1).style = 'popupmenu';
    Formats(end,1).format = 'input';
    Formats(end,1).items = {'FASTER','TAPEEG','disabled'};
    Formats(end,1).callback = {@set_AVGroutine,'@ALL','@ALL',header,cfg};
    
    Prompt(end+1,:) = {'settings','BAD',''};
    Formats(end+1,1).type = 'check';
    Formats(end,1).callback = {@set_AVGroutine,'@ALL','@ALL',header,cfg};
    
    Formats(end+1,1).type = 'none';
    Formats(end,1).span = [1 3];
    
    Prompt(end+1,:) = {'Method','method',''};
    Formats(end+1,1).type = 'list';
    Formats(end,1).style = 'popupmenu';
    Formats(end,1).format = 'input';
    Formats(end,1).items = {'single','percent'};
    Formats(end,1).callback = {@set_AVGmethod,'percent','percent','method'};
    
    Prompt(end+1,:) = {'','percent','% (channels)'};
    Formats(end+1,1).type = 'edit';
    Formats(end,1).format = 'float';
    Formats(end,1).limits = [0 inf];
    Formats(end,1).size = 25;
    
    Prompt(end+1,:) = {' ','','%'};
    Formats(end+1,1).type = 'text';
    
    [settings,Cancelled] = inputsdlg(Prompt,'Reject',Formats,settings);
    if Cancelled == true
        settings = [];
    end
end
   
function settings = set_AVGroutine(settings,header,cfg)
    if strcmp(settings.routine,'FASTER')
        if ~isempty(settings.BAD)
            settings.BAD = 1.5;
        end
        Prompt = cell(0,3);
        Formats = [];
        Prompt(end+1,:) = {'Standard deviation for good/bad sweep','BAD',''};
        Formats(end+1,1).type = 'edit';
        Formats(end,1).format = 'float';
        Formats(end,1).limits = [0 inf];
        Formats(end,1).size = 25;
        Formats(end,1).span = [1 3];
        [settings,Cancelled] = inputsdlg(Prompt,'Reject',Formats,settings);
        if Cancelled == 1
            settings.BAD = [];
            settings.routine = 'disabled';
        end
    elseif strcmp(settings.routine,'TAPEEG')
        if isnumeric(settings.BAD) | ~isfield(settings.BAD,'percentbad')
            settings.BAD.percentbad = 1;
            settings.BAD.freqlim50 = [];
            settings.BAD.freqlim60 = [];
            settings.BAD.freqlimlow = [];
            settings.BAD.freqlimhigh = [];
            settings.BAD.spectslow = [0.5 2];
            settings.BAD.spectshigh = [15 50];
            settings.BAD.zvaluebroken = [];
            settings.BAD.zvaluevars = 3.5;
            settings.BAD.zvaluehurst = 3.5;
            settings.BAD.zvaluemedian = 3.5;
            settings.BAD.zvalueamplitude = 3.5;
            if ~isfield(header,'numdatachannels') | header.numdatachannels > 72
                settings.BAD.zvaluecorr = 3.5;
            else
                settings.BAD.zvaluecorr = 0;
            end
            settings.BAD.LAPL.lap_maxdistance = 2.5;
            settings.BAD.LAPL.lap_weightmaxdistance = 30;
            settings.BAD.LAPL.lap_excluderef = true;
        end
        settings.BAD = lab_set_detect_bad(settings.BAD,header,cfg,0,0,0,0,0,0,0,0,1,0);
        if isempty(settings.BAD)
            settings.routine = 'disabled';
        end
    else
        settings.BAD = [];
        settings.routine = 'disabled';
    end    
end

function percent = set_AVGmethod(percent,method)
    if strcmp(method,'single')
        percent = [];
    elseif strcmp(method,'percent') & isempty(percent)
        percent = 5;
    end
end

function Correct = set_correct(Correct,header,cfg)
    if ~isfield(Correct,'method')
        Correct.method = 'spherical';
    end
    if ~isfield(Correct,'BAD') | isempty(Correct.BAD)
        Correct.BAD.freqlim50 = 0;
        Correct.BAD.freqlim60 = 0;
        Correct.BAD.freqlimlow = 0;
        Correct.BAD.freqlimhigh = 0;
        Correct.BAD.spectslow = [0.5 2];
        Correct.BAD.spectshigh = [15 50];
        Correct.BAD.zvaluevars = 2.5;
        Correct.BAD.zvaluehurst = 2.5;
        Correct.BAD.zvaluemedian = 0;
        Correct.BAD.zvaluecorr = 2.5;
        Correct.BAD.LAPL.lap_maxdistance = 2.5;
        Correct.BAD.LAPL.lap_weightmaxdistance = 30;
        Correct.BAD.LAPL.lap_excluderef = true;
    end
    
    Prompt = cell(0,2);
    Formats = [];
    
    Prompt(end+1,:) = {'Detect bad','BAD'};
    Formats(end+1,1).type = 'edit';
    Formats(end,1).format = 'result';
    Formats(end,1).callback = {@lab_set_detect_bad,'BAD','BAD',header,cfg,0,0,0,0,0,0,0,0,0,0};
    Formats(end,1).size = [150 160];
    
    Formats(end+1,1).type = 'none';
    
    Prompt(end+1,:) = {'Interpolate method','method'};
    Formats(end+1,1).type = 'list';
    Formats(end,1).style = 'popupmenu';
    Formats(end,1).format = 'input';
    Formats(end,1).items = {'disabled';'spherical';'3D'};
    
    [Correct,Cancelled] = inputsdlg(Prompt,'Detect bad/Interpolate',Formats,Correct);
    if Cancelled == 1
        Correct = [];
    end
end

function correctbaseline = set_correctbaseline(correctbaseline,markerOffset)
    if strcmp(correctbaseline,'pre-marker period') & markerOffset <= 0
        correctbaseline = 'all TFs';
    end
end

function correctbaseline = set_markerOffset(correctbaseline,markerOffset)
    if markerOffset <= 0
        correctbaseline = false;
    end
end