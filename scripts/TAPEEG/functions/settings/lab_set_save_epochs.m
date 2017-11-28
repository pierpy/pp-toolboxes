function [cfg,skipprocessing] = lab_set_save_epochs(cfg,header,doshort,domarker,dofreqs,domicrocorr)

skipprocessing = 0;

if ~exist('domicrocorr','var')
    domicrocorr = false;
end
if ~exist('dofreqs','var')
    dofreqs = true;
end
if ~exist('domarker','var')
    domarker = true;
end
if ~exist('doshort','var')
    doshort = false;
end

global THeader
if ~isempty(THeader)
    header = THeader;
elseif ~exist('header','var')
    header = [];
end

if ~exist('cfg','var') | ~isfield(cfg,'EPOCH') | ~isfield(cfg.EPOCH,'minimalpart')
    cfg.EPOCH.folder = 'Epochs';
    cfg.EPOCH.minimalpart = 12;
    cfg.EPOCH.percentgood = 100;
    cfg.EPOCH.markerexclude = {''};
    cfg.EPOCH.markerinclude = {''};
    cfg.EPOCH.markerstart ='';
    cfg.EPOCH.markerstop = '';
    cfg.EPOCH.replacebad = false;
    cfg.EPOCH.interpolate3D = false;
    if dofreqs == true
        cfg.EPOCH.BAD.freqlim50 = 50;
        cfg.EPOCH.BAD.freqlim60 = 50;
        cfg.EPOCH.BAD.freqlimlow = 70;
        cfg.EPOCH.BAD.freqlimhigh = 50;
    else
        cfg.EPOCH.BAD.freqlim50 = [];
        cfg.EPOCH.BAD.freqlim60 = [];
        cfg.EPOCH.BAD.freqlimlow = [];
        cfg.EPOCH.BAD.freqlimhigh = [];
    end
    cfg.EPOCH.BAD.spectslow = [0.5 2];
    cfg.EPOCH.BAD.spectshigh = [15 50];
    cfg.EPOCH.BAD.zvaluebroken = 0;
    cfg.EPOCH.BAD.zvaluevars = 4;
    cfg.EPOCH.BAD.zvaluehurst = 3;
    cfg.EPOCH.BAD.zvaluemedian = 3;
    if ~isfield(header,'numdatachannels') | header.numdatachannels > 72
        cfg.EPOCH.BAD.zvaluecorr = 3;
    else
        cfg.EPOCH.BAD.zvaluecorr = 0;
    end
    cfg.EPOCH.BAD.LAPL.lap_maxdistance = 2.5;
    cfg.EPOCH.BAD.LAPL.lap_weightmaxdistance = 30;
    cfg.EPOCH.BAD.LAPL.lap_excluderef = true;
    cfg.EPOCH.BAD.PEAK2MIN.lowfreqpeak = 4;
    cfg.EPOCH.BAD.PEAK2MIN.highfreqpeak = 14;
    cfg.EPOCH.BAD.PEAK2MIN.MinPeak2Min = 1.3;
    cfg.EPOCH.BAD.PEAK2MIN.mode = 'weighted';
    cfg.EPOCH.BAD.PEAK2MIN.threshold = [];
    cfg.EPOCH.BAD.PEAK2MIN.factor = 3;
    if isfield(header,'numdatachannels')
        cfg.EPOCH.BAD.PEAK2MIN.BAchannels = lab_get_bactivity(header.numdatachannels);
    elseif isfield(header,'numchannels')
        cfg.EPOCH.BAD.PEAK2MIN.BAchannels = lab_get_bactivity(header.numchannels);
    else
        cfg.EPOCH.BAD.PEAK2MIN.BAchannels = [];
    end
end
if domicrocorr == true & ~isfield(cfg.EPOCH.BAD,'MicroCorr')
    cfg.EPOCH.BAD.MicroCorr = true;
elseif isfield(cfg.EPOCH.BAD,'MicroCorr')
    cfg.EPOCH.BAD = rmfield(cfg.EPOCH.BAD,'MicroCorr');
end
if doshort == false & ~isfield(cfg.EPOCH,'eegsource')
    cfg.EPOCH.eegsource = 'input';
    cfg.EPOCH.REF = [];
    cfg.EPOCH.montage = [];
    cfg.EPOCH.length = 4096;
    cfg.EPOCH.eformat = {'eph'};
    cfg.EPOCH.interpolatebad = true;
    cfg.EPOCH.plotquality = true;
end

Marker = lab_getall_markers(header,cfg);

Prompt = cell(0,2);
Formats = {};
if doshort == false
    Prompt(end+1,:) = {'Output-folder', 'folder'};
    Formats(end+1,1).type = 'edit';
    Formats(end,1).format = 'text';
    Formats(end,1).size = 100;
    Formats(end,1).span = [1 3];
    
    Formats(end+1,1).type = 'none';
    Formats(end,1).span = [1 3];
    
    Prompt(end+1,:) = {'Reference','eegsource'};
    Formats(end+1,1).type = 'list';
    Formats(end,1).style = 'popupmenu';
    Formats(end,1).format = 'input';
    if ~max(strcmp({'channels','mean','median','laplacian','montage','input'},cfg.EPOCH.eegsource))
        Formats(end,1).items = {cfg.EPOCH.eegsource,'channels','mean','median','laplacian','montage','input'};
    else
        Formats(end,1).items = {'channels','mean','median','laplacian','montage','input'};
    end
    Formats(end,1).callback = {@lab_get_eegsource,'@ALL','@ALL',cfg,header};
    
    Prompt(end+1,:) = {'Montage','montage'};
    Formats(end+1,1).type = 'check';
    Formats(end,1).callback = {@lab_load_montage,'montage','montage',cfg,header};
    
    Prompt(end+1,:) = {'Laplacian','LAPL'};
    Formats(end+1,1).type = 'check';
    Formats(end,1).callback = {@lab_get_laplacian,'LAPL','LAPL'};
    
    Prompt(end+1,:) = {'Interpolate bad channels','interpolatebad'};
    Formats(end+1,1).type = 'check';
    Formats(end,1).size = [0 -1];
    Formats(end,1).span = [1 3];
    
    Formats(end+1,1).type = 'none';
    Formats(end,1).span = [1 3];
    
    Prompt(end+1,:) = {'Epoch length (samples)','length'};
    Formats(end+1,1).type = 'edit';
    Formats(end,1).format = 'integer';
    Formats(end,1).limits = [0 99999];
    Formats(end,1).size = 50;
    Formats(end,1).span = [1 3];
end

Prompt(end+1,:) = {['Minimal part of summed epochs on eeg length' sprintf('\n')  ...
    '(between 0 and 1 / >1 for fixed number of epochs)'],'minimalpart'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'float';
Formats(end,1).limits = [0 inf];
Formats(end,1).size = 35;
Formats(end,1).span = [1 3];

Formats(end+1,1).type = 'none';
Formats(end,1).span = [1 3];

Prompt(end+1,:) = {'Detect bad','BAD'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'result';
Formats(end,1).callback = {@lab_set_detect_bad,'BAD','BAD',header,cfg,0,0,0,0,0,1,0,0,0,0};
Formats(end,1).size = [150 200];
Formats(end,1).span = [1 3];

Prompt(end+1,:) = {'Minimal percent of good channels in single epoch','percentgood'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'float';
Formats(end,1).limits = [0 100];
Formats(end,1).size = 35;
Formats(end,1).span = [1 3];

Formats(end+1,1).type = 'none';
Formats(end,1).span = [1 3];

if domarker == true
    Prompt(end+1,:) = {'Markers to exclude (strings / * / all)', 'markerexclude'};
    Formats(end+1,1).type = 'edit';
    Formats(end,1).format = 'vector';
    Formats(end,1).size = 220;
    Formats(end,1).span = [1 3];
    if ~isempty(Marker)
        Formats(end,1).items = [{'all'} {'*'} Marker(:)'];
    end
    
    Prompt(end+1,:) = {'Markers to include  (strings / * / all)', 'markerinclude'};
    Formats(end+1,1).type = 'edit';
    Formats(end,1).format = 'vector';
    Formats(end,1).size = 220;
    Formats(end,1).span = [1 3];
    if ~isempty(Marker)
        Formats(end,1).items = [{'all'} {'*'} Marker(:)'];
    end
end

Prompt(end+1,:) = {'Marker to start taking epochs (empty / string / *)', 'markerstart'};
Formats(end+1,1).size = 100;
Formats(end,1).span = [1 3];
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
Formats(end,1).span = [1 3];
if ~isempty(Marker)
    Formats(end,1).type = 'list';
    Formats(end,1).style = 'popupmenu';
    Formats(end,1).format = 'input';
    Formats(end,1).items = [{''} {'all'} {'*'} Marker(:)'];
else
    Formats(end,1).type = 'edit';
    Formats(end,1).format = 'text';
end

if doshort == false
    Prompt(end+1,:) = {' ',''};
    Formats(end+1,1).type = 'text';
    Formats(end,1).span = [1 3];
    
    Prompt(end+1,:) = {'File format','eformat'};
    Formats(end+1,1).type = 'list';
    Formats(end,1).style = 'listbox';
    Formats(end,1).format = 'input';
    Formats(end,1).items = {'ris';'sef';'edf';'eph';'ep';'txt'};
    Formats(end,1).limits = [0 2]; % multi-select
    Formats(end,1).size = [60 100];
    Formats(end,1).span = [4 1];
    Formats(end,1).callback = {@lab_get_format,'scaletxt','eformat','scaletxt'};
    
    Prompt(end+1,:) = {'Replace bad channels','replacebad'};
    Formats(end+1,1).type = 'check';
    Formats(end,1).size = [-1 -1];
    Formats(end,1).span = [1 2];
    
    Prompt(end+1,:) = {'3D-Interpolation of bad channels','interpolate3D'};
    Formats(end+1,1).type = 'check';
    Formats(end,1).size = [-1 -1];
    Formats(end,1).span = [1 2];

    Prompt(end+1,:) = {'Plot quality control','plotquality'};
    Formats(end+1,1).type = 'check';
    Formats(end,1).span = [1 2];
    Formats(end,1).size = [-1 -1];
    
    Prompt(end+1,:) = {'Scale TXT','scaletxt'};
    Formats(end+1,1).type = 'check';
    Formats(end,1).span = [1 2];
    Formats(end,1).callback = {@lab_show_result,'scaletxt','scaletxt'};
end

[cfg.EPOCH,Cancelled] = inputsdlg(Prompt,'Epoch settings',Formats,cfg.EPOCH,2);
if isempty(cfg.EPOCH) | Cancelled == 1
    cfg.EPOCH= [];
    skipprocessing = 1;
end
