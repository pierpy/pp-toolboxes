function [settings,skipprocessing] = lab_set_segments_auto(settings,header,cfg)

skipprocessing = 0;

if ~exist('cfg','var')
    cfg = [];
end
if ~exist('header','var')
    header = [];
end

if ~exist('settings','var')  | isempty(settings)
    settings.mergefolder = false;
    settings.minduration = 35;
    settings.durationall = 180;
    settings.percentgood = 100;
    settings.markerexclude = {''};
    settings.markerinclude = {''};
    settings.markerstart ='';
    settings.markerstop = '';
    settings.interpolate3D = false;
    settings.BAD.length = 2;
    settings.BAD.freqlim50 = 50;
    settings.BAD.freqlim60 = 50;
    settings.BAD.freqlimlow = 70;
    settings.BAD.freqlimhigh = 50;
    settings.BAD.spectslow = [0.5 2];
    settings.BAD.spectshigh = [15 50];
    settings.BAD.zvaluebroken = 4;
    settings.BAD.zvaluevars = 3;
    settings.BAD.zvaluehurst = 3;
    settings.BAD.zvaluemedian = 0;
    if ~isfield(header,'numdatachannels') | header.numdatachannels > 72
        settings.BAD.zvaluecorr = 3;
    else
        settings.BAD.zvaluecorr = 0;
    end
    settings.BAD.LAPL.lap_maxdistance = 2.5;
    settings.BAD.LAPL.lap_weightmaxdistance = 30;
    settings.BAD.LAPL.lap_excluderef = true;
    settings.BAD.PEAK2MIN.lowfreqpeak = 4;
    settings.BAD.PEAK2MIN.highfreqpeak = 14;
    settings.BAD.PEAK2MIN.MinPeak2Min = 1.3;
    settings.BAD.PEAK2MIN.mode = 'threshold';
    settings.BAD.PEAK2MIN.factor = 2;
    settings.BAD.PEAK2MIN.threshold = 1.5;
    if isfield(header,'numdatachannels')
        settings.BAD.PEAK2MIN.BAchannels = lab_get_bactivity(header.numdatachannels);
    elseif isfield(header,'numchannels')
        settings.BAD.PEAK2MIN.BAchannels = lab_get_bactivity(header.numchannels);
    else
        settings.BAD.PEAK2MIN.BAchannels = [];
    end
end

Marker = lab_getall_markers(header,cfg);

Prompt = cell(0,2);
Formats = {};
    
Prompt(end+1,:) = {'Minimal duration of single segment (sec)','minduration'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'integer';
Formats(end,1).limits = [0 99999];
Formats(end,1).size = 50;

Prompt(end+1,:) = {'Merge EEG''s in same folder','mergefolder'};
Formats(end+1,1).type = 'check';
    
Prompt(end+1,:) = {'Duration of all segments (sec)','durationall'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'integer';
Formats(end,1).limits = [0 99999];
Formats(end,1).size = 50;
Formats(end,1).span = [1 2];

Formats(end+1,1).type = 'none';
Formats(end+1,1).type = 'none';

Prompt(end+1,:) = {'Percent of good channels in single segment','percentgood'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'float';
Formats(end,1).limits = [0 100];
Formats(end,1).size = 35;
Formats(end,1).span = [1 2];

Prompt(end+1,:) = {'Detect segments','BAD'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'result';
Formats(end,1).callback = {@lab_set_detect_bad,'BAD','BAD',header,cfg,0,0,1,0,0,1,0,0,0,0};
Formats(end,1).size = [150 200];
Formats(end,1).span = [1 2];

Formats(end+1,1).type = 'none';
Formats(end+1,1).type = 'none';

Prompt(end+1,:) = {'Markers to exclude (strings / * / all)', 'markerexclude'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'vector';
Formats(end,1).size = 220;
Formats(end,1).span = [1 2];
if ~isempty(Marker)
    Formats(end,1).items = [{'all'} {'*'} Marker(:)'];
end

Prompt(end+1,:) = {'Markers to include  (strings / * / all)', 'markerinclude'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'vector';
Formats(end,1).size = 220;
Formats(end,1).span = [1 2];
if ~isempty(Marker)
    Formats(end,1).items = [{'all'} {'*'} Marker(:)'];
end

Prompt(end+1,:) = {'Marker to start taking epochs (empty / string / *)', 'markerstart'};
Formats(end+1,1).size = 100;
Formats(end,1).span = [1 2];
if ~isempty(Marker)
    Formats(end,1).type = 'list';
    Formats(end,1).style = 'popupmenu';
    Formats(end,1).format = 'input';
    Formats(end,1).items = [{''} {'all'} {'*'} Marker(:)'];
else
    Formats(end,1).type = 'edit';
    Formats(end,1).format = 'text';
end

Prompt(end+1,:) = {'Marker to stop taking epochs (empty / string / *)', 'markerstop'};
Formats(end+1,1).size = 100;
Formats(end,1).span = [1 2];
if ~isempty(Marker)
    Formats(end,1).type = 'list';
    Formats(end,1).style = 'popupmenu';
    Formats(end,1).format = 'input';
    Formats(end,1).items = [{''} {'all'} {'*'} Marker(:)'];
else
    Formats(end,1).type = 'edit';
    Formats(end,1).format = 'text';
end

Prompt(end+1,:) = {'Montage for export','montage'};
Formats(end+1,1).type = 'check';
Formats(end,1).callback = {@lab_load_montage,'montage','montage',cfg,header,'MontageEDF.xls'};

[settings,Cancelled] = inputsdlg(Prompt,'Auto Segments',Formats,settings,2);
if isempty(settings) | Cancelled == 1
    settings = [];
    skipprocessing = 1;
    return
end

