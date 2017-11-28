% settings for detect bad
%
% [settings,skipprocessing] = lab_set_detect_bad(settings,header,cfg,doeye,doecg,doepoch,doica,domark,dopeak,doavg,dofile,domarker)
%
% F. Hatz 2013

function [settings,skipprocessing] = lab_set_detect_bad(settings,header,cfg,doeye,doecg,doperiod,doica,domark,dopeak,doavg,dofile,domarker,dofixed)

skipprocessing = 0;
if ~exist('dofixed','var')
    dofixed = 1;
end
if ~exist('domarker','var')
    domarker = 1;
end
if ~exist('dofile','var')
    dofile = 1;
end
if ~exist('doavg','var')
    doavg = 1;
end
if ~exist('dopeak','var')
    dopeak = 1;
end
if ~exist('domark','var')
    domark = 1;
end
if ~exist('doica','var')
    doica = 2;
end
if ~exist('doperiod','var')
    doperiod = 2;
end
if ~exist('doecg','var')
    doecg = 1;
end
if ~exist('doeye','var')
    doeye = 1;
end
if ~exist('cfg','var')
    cfg = [];
end
if ~exist('header','var')
    header = [];
end
if ~exist('settings','var')
    settings = [];
end

Marker = lab_getall_markers(header,cfg);

if isempty(settings)
    if dofile == 1
        settings.filemethod = 'add';
    end
    if doperiod > 0
        settings.length = 4;
    end
    if doperiod == 2
        settings.percentbad = 70;
    end
    settings.freqlim50 = 50;
    settings.freqlim60 = 50;
    settings.freqlimlow = 70;
    settings.freqlimhigh = 50;
    settings.spectslow = [0.5 2];
    settings.spectshigh = [15 50];
    if doica == 1 | doperiod == 0
        settings.zvaluebroken = 0;
    else
        settings.zvaluebroken = 4;
    end
    if doica == 1
        settings.zvaluevars = 0;
    elseif isfield(header,'interpolated') & ~isempty(header.interpolated)
        settings.zvaluevars = 4;
    else
        settings.zvaluevars = 3;
    end
    settings.zvaluehurst = 3;
    if doica == 1
        settings.zvaluemedian = 3;
        settings.zvaluekurtosis = 3;
    else
        settings.zvaluemedian = 0;
        if ~isfield(header,'numdatachannels') | header.numdatachannels > 72
            settings.zvaluecorr = 3;
        else
            settings.zvaluecorr = 0;
        end
        settings.LAPL.lap_maxdistance = 2.5;
        settings.LAPL.lap_weightmaxdistance = 30;
        settings.LAPL.lap_excluderef = true;
    end
    settings.amplitude = 0;
    if dopeak == 1
        if doica == 0
            settings.PEAK2MIN.lowfreqpeak = 4;
            settings.PEAK2MIN.highfreqpeak = 14;
            settings.PEAK2MIN.MinPeak2Min = 1.3;
            settings.PEAK2MIN.mode = 'threshold';
            settings.PEAK2MIN.factor = 2;
            settings.PEAK2MIN.threshold = 1.5;
        else
            settings.PEAK2MIN = [];
        end
    end
    if doavg == 1
        settings.AVG = [];
        settings.AVGstd = 2.5;
        settings.AVGmode = 'detect bad';
    end
end
if isfield(settings,'MicroCorr')
    domicrocorr = true;
else
    domicrocorr = false;
end

if doeye > 0 & ~isfield(settings,'eog')
    if exist('header','var') & isfield(header,'eog_ch') & header.eog_ch > 0
        settings.zvalueeye = 3;
        settings.eog = header.eog_ch(:);
    elseif isfield(cfg,'EXTRA') & isfield(cfg.EXTRA,'numdatachans') & cfg.EXTRA.numdatachans == 257
        settings.zvalueeye = 3;
        settings.eog = [37,241,46,244;18,238,10,234;31,252,32,253;31,226,25,225];
    elseif isfield(header,'numdatachannels') & header.numdatachannels == 257
        settings.zvalueeye = 3;
        settings.eog = [37,241,46,244;18,238,10,234;31,252,32,253;31,226,25,225];
    else
        settings.zvalueeye = 3;
        settings.eog = [];
    end
end
if doecg > 0 & ~isfield(settings,'ecgdetect')
    if exist('header','var') & isfield(header,'ecg_ch') & ~isempty(header.ecg_ch) & header.ecg_ch ~= 0
        settings.ecgdetect = 1;
        settings.ecg_ch = header.ecg_ch;
    else
        settings.ecgdetect = 2;
        settings.ecg_ch = [50 120];
    end
end
if domark > 0 & ~isfield(settings,'markbad')
    settings.markbad = false;
    settings.markbadValue = 25;
end
if dofixed > 0 & ~isfield(settings,'fixedbad')
    settings.fixedbad = [];
end

if isfield(cfg,'MAIN') & isfield(cfg.MAIN,'auto') & cfg.MAIN.auto == 1
    return
end

Prompt = cell(0,2);
Formats = {};

if dofile > 0
    Prompt(end+1,:) = {'Read bad channels from file (.txt,.xls,.info)','filemethod'};
    Formats(end+1,1).type = 'list';
    Formats(end,1).style = 'popupmenu';
    Formats(end,1).format = 'input';
    Formats(end,1).items = {'add','only','skip','or'};
    Formats(end,1).span = [1 3];
end
if dofixed > 0
    Prompt(end+1,:) = {'Set fixed bad channels','fixedbad'};
    Formats(end+1,1).type = 'edit';
    Formats(end,1).format = 'vector';
    Formats(end,1).limits = [-inf inf];
    Formats(end,1).size = 120;
    Formats(end,1).span = [1 3];
end

if ~isempty(Prompt)
    Formats(end+1,1).type = 'none';
    Formats(end,1).span = [1 3];
end

if doperiod > 0
    Prompt(end+1,:) = {'Period length for analysis (seconds)','length'};
    Formats(end+1,1).type = 'edit';
    Formats(end,1).format = 'float';
    Formats(end,1).limits = [0 9999];
    Formats(end,1).size = 30;
end
if doperiod == 2
    Prompt(end+1,:) = {'Allowed percent bad periods per channel','percentbad'};
    Formats(end+1,1).type = 'edit';
    Formats(end,1).format = 'float';
    Formats(end,1).limits = [0 100];
    Formats(end,1).size = 30;
    Formats(end,1).span = [1 2];
end

if ~isempty(Prompt)
    Formats(end+1,1).type = 'none';
    Formats(end,1).span = [1 3];
end

Prompt(end+1,:) = {'Frequencies (max percent):',''};
Formats(end+1,1).type = 'text';
Formats(end,1).span = [1 3];

Prompt(end+1,:) = {'50Hz (0=off)','freqlim50'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'float';
Formats(end,1).limits = [0 100];
Formats(end,1).size = 40;
Formats(end,1).span = [1 3];

Prompt(end+1,:) = {'60Hz (0=off)','freqlim60'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'float';
Formats(end,1).limits = [0 100];
Formats(end,1).size = 40;
Formats(end,1).span = [1 3];

Prompt(end+1,:) = {'spect-range low (0=off)','freqlimlow'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'float';
Formats(end,1).limits = [0 100];
Formats(end,1).size = 40;

Prompt(end+1,:) = {'Range (lowfreq highfreq)','spectslow'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'vector';
Formats(end,1).limits = [-inf inf];
Formats(end,1).size = 50;
Formats(end,1).span = [1 2];

Prompt(end+1,:) = {'spect-range high (0=off)','freqlimhigh'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'float';
Formats(end,1).limits = [0 100];
Formats(end,1).size = 40;

Prompt(end+1,:) = {'Range (lowfreq highfreq)','spectshigh'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'vector';
Formats(end,1).limits = [-inf inf];
Formats(end,1).size = 50;
Formats(end,1).span = [1 2];

Formats(end+1,1).type = 'none';
Formats(end,1).span = [1 3];

Prompt(end+1,:) = {'Z-value broken (0=off)','zvaluebroken'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'float';
Formats(end,1).limits = [0 10];
Formats(end,1).size = 40;
Formats(end,1).span = [1 3];

Prompt(end+1,:) = {'Z-value vars (0=off)','zvaluevars'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'float';
Formats(end,1).limits = [0 10];
Formats(end,1).size = 40;
Formats(end,1).span = [1 3];

Prompt(end+1,:) = {'Z-value hurst (0=off)','zvaluehurst'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'float';
Formats(end,1).limits = [0 10];
Formats(end,1).size = 40;
Formats(end,1).span = [1 3];

Prompt(end+1,:) = {'Z-value median (0=off)','zvaluemedian'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'float';
Formats(end,1).limits = [0 10];
Formats(end,1).size = 40;
Formats(end,1).span = [1 3];

Prompt(end+1,:) = {'Z-value amplitude (0=off)','zvalueamplitude'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'float';
Formats(end,1).limits = [0 10];
Formats(end,1).size = 40;
Formats(end,1).span = [1 3];

if doica > 0
    Prompt(end+1,:) = {'Z-value kurtosis (0=off)','zvaluekurtosis'};
    Formats(end+1,1).type = 'edit';
    Formats(end,1).format = 'float';
    Formats(end,1).limits = [0 10];
    Formats(end,1).size = 40;
    Formats(end,1).span = [1 3];
end

if doica == 0 | doica == 2
    Prompt(end+1,:) = {'Z-value topo-corr (0=off)','zvaluecorr'};
    Formats(end+1,1).type = 'edit';
    Formats(end,1).format = 'float';
    Formats(end,1).limits = [0 10];
    Formats(end,1).size = 40;
    Formats(end,1).callback = {@set_corr,'@ALL','@ALL'};
    
    Prompt(end+1,:) = {'Laplacian','LAPL'};
    Formats(end+1,1).type = 'check';
    Formats(end,1).callback = {@lab_get_laplacian,'LAPL','LAPL',1};
    Formats(end,1).span = [1 2];
end

if doeye > 0
    Formats(end+1,1).type = 'none';
    Formats(end,1).span = [1 3];
    
    Prompt(end+1,:) = {'Z-value eye (0=off)','zvalueeye'};
    Formats(end+1,1).type = 'edit';
    Formats(end,1).format = 'float';
    Formats(end,1).limits = [0 10];
    Formats(end,1).size = 40;
    Formats(end,1).callback = {@correct_eog,'eog','eog','zvalueeye',cfg,header};
    
    Prompt(end+1,:) = {'Channel-list','eog'};
    Formats(end+1,1).type = 'edit';
    Formats(end,1).format = 'result';
    Formats(end,1).callback = {@lab_get_eog,'eog','eog',header,cfg};
    Formats(end,1).size = 70;
    Formats(end,1).span = [1 2];
end

if doecg > 0
    Formats(end+1,1).type = 'none';
    Formats(end,1).span = [1 3];

    Prompt(end+1,:) = {'Ecg artifact','ecgdetect'};
    Formats(end+1,1).type = 'list';
    Formats(end,1).style = 'popupmenu';
    Formats(end,1).items = {'ecg_channel','heartbeat','disabled'};
    Formats(end,1).callback = {@correct_ecg,'ecg_ch','ecg_ch','ecgdetect',cfg,header};
    
    Prompt(end+1,:) = {'Ecg channel / heartbeat range','ecg_ch'};
    Formats(end+1,1).type = 'edit';
    Formats(end,1).format = 'vector';
    Formats(end,1).limits = [0 inf];
    Formats(end,1).size = 60;
    Formats(end,1).span = [1 2];
end

if dopeak > 0
    Formats(end+1,1).type = 'none';
    Formats(end,1).span = [1 3];
    
    Prompt(end+1,:) = {'Peak2Min','PEAK2MIN'};
    Formats(end+1,1).type = 'check';
    Formats(end,1).callback = {@lab_get_peak2min,'PEAK2MIN','PEAK2MIN',header,cfg};
    Formats(end,1).span = [1 3];
end

if domarker > 0
    Formats(end+1,1).type = 'none';
    Formats(end,1).span = [1 3];
    
    Prompt(end+1,:) = {'Markers to exclude (strings / * / all)',''};
    Formats(end+1,1).type = 'text';
    Formats(end,1).span = [1 3];
    
    Prompt(end+1,:) = {'','markerexclude'};
    Formats(end+1,1).type = 'edit';
    Formats(end,1).format = 'vector';
    Formats(end,1).limits = [0 1];
    Formats(end,1).size = 200;
    if ~isempty(Marker)
        Formats(end,1).items = [{'all'} Marker(:)'];
        Formats(end,1).callback = {@lab_get_marker,'markerexclude','markerexclude'};
    end
    Formats(end,1).span = [1 3];
end

if doavg > 0
    Formats(end+1,1).type = 'none';
    Formats(end,1).span = [1 3];
    
    Prompt(end+1,:) = {'Standard deviation in average','AVG'};
    Formats(end+1,1).type = 'check';
    Formats(end,1).callback = {@set_average,'@ALL','@ALL',header};
    
    Prompt(end+1,:) = {'Std >','AVGstd'};
    Formats(end+1,1).type = 'edit';
    Formats(end,1).format = 'float';
    Formats(end,1).limits = [0 25];
    Formats(end,1).size = 40;
    
    Prompt(end+1,:) = {'','AVGmode'};
    Formats(end+1,1).type = 'list';
    Formats(end,1).style = 'popupmenu';
    Formats(end,1).format = 'input';
    Formats(end,1).items = {'detect bad','correct bad'};
end

if domark > 0
    Formats(end+1,1).type = 'none';
    Formats(end,1).span = [1 3];

    Prompt(end+1,:) = {'Create markers for bad periods' 'markbad'};
    Formats(end+1,1).type = 'check';
    Formats(end,1).format = 'integer';
    Formats(end,1).callback = {@set_mark,'@ALL','@ALL'};
    
    Prompt(end+1,:) = {'Percent of bad channels','markbadValue'};
    Formats(end+1,1).type = 'edit';
    Formats(end,1).format = 'float';
    Formats(end,1).limits = [0 100];
    Formats(end,1).size = 40;
    Formats(end,1).span = [1 2];
end

if domicrocorr == true
    Prompt(end+1,:) = {'Microstates-Correlations' 'MicroCorr'};
    Formats(end+1,1).type = 'check';
    
    Prompt(end+1,:) = {'Take negativ values' 'MicroCorrNeg'};
    Formats(end+1,1).type = 'check';
    Formats(end,1).span = [1 2];
end

[settings,Cancelled] = inputsdlg(Prompt,'Detect Bad',Formats,settings);
if isempty(settings) | Cancelled == 1
    settings = [];
    skipprocessing = 1;
    return
end
pause(0.2);

end

function ecg_ch = correct_ecg(ecg_ch,ecgdetect,cfg,header)
    if ecgdetect == 1
        if isfield(cfg,'EXTRA') & isfield(cfg.EXTRA,'ecg_ch')
            ecg_ch = cfg.EXTRA.ecg_ch;
        elseif isfield(header,'ecg_ch')
            ecg_ch = header.ecg_ch;
        elseif length(ecg_ch) > 1 | (isfield(header,'numchannels') & ecg_ch > header.numchannels)
            ecg_ch = [];
        end
    elseif ecgdetect == 2 & length(ecg_ch) ~= 2
        ecg_ch = [50 120];
    elseif ecgdetect == 3
        ecg_ch = [];
    end
end

function eog = correct_eog(eog,zvalueeye,cfg,header)
    if zvalueeye > 0 & isempty(eog)
        if exist('header','var') & isfield(header,'eog_ch') & header.eog_ch > 0
            eog = header.eog_ch;
        elseif isfield(cfg,'EXTRA') & isfield(cfg.EXTRA,'numdatachans') & cfg.EXTRA.numdatachans == 257
            eog = [37,241,46,244;18,238,10,234;31,252,32,253;31,226,25,225];
        elseif isfield(header,'numdatachannels') & header.numdatachannels == 257
            eog = [37,241,46,244;18,238,10,234;31,252,32,253;31,226,25,225];
        end
    end
end

function settings = set_average(settings,header)
    if isempty(settings.AVG) | ~isfield(settings,'AVG')
        settings.AVG.marker = {'*'};
        settings.AVG.combinemarker = true;
        settings.AVG.markerlength = 500;
        settings.AVG.markerOffset = 100;
        settings.AVG.AVGmethod = 'mean';
        settings.AVG.Reject.routine = 'FASTER';
        settings.AVG.Reject.BAD = 1.5;
        settings.AVG.Reject.method = 'single';
        settings.AVG.Reject.percent = [];
    end
    settings = lab_set_average(settings,header,[],1);
    if ~isempty(settings.AVG) & isempty(settings.AVGstd)
        settings.AVGstd = 2.5;
    elseif isempty(settings.AVG)
        settings.AVGstd = [];
    end
end

function settings = set_mark(settings)
    if settings.markbad == false
        settings.markbad = true;
        settings.length = 1;
    else
        settings.markbad = false;
        if settings.length == 1
            settings.length = 4;
        end
    end
end

function settings = set_corr(settings)
    if settings.zvaluecorr > 0 & isempty(settings.LAPL)
        settings.LAPL.lap_maxdistance = 2.5;
        settings.LAPL.lap_weightmaxdistance = 30;
        settings.LAPL.lap_excluderef = true;
    elseif settings.zvaluecorr == 0 | isempty(settings.zvaluecorr)
        settings.LAPL = [];
    end
end