function [cfg,skipprocessing] = lab_set_remove_dbs(cfg,header)

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

if ~isfield(cfg,'STIM') | isempty(cfg.STIM)
    cfg.STIM.folder = 'RemoveDBS';
    cfg.STIM.DETECT.stimrange = [100,150;80,100;150,250];
    cfg.STIM.DETECT.stimmaxpercent = 70;
    cfg.STIM.PCA = true;
    cfg.STIM.ICA.type = 'runica';
    cfg.STIM.ICA.extended = 0;
    cfg.STIM.ICA.pca = 0;
    cfg.STIM.ICA.good = true;
    cfg.STIM.ICA.epoch = 0;
    cfg.STIM.ICA.plottopo = false;
    cfg.STIM.ICA.BAD = [];
    cfg.STIM.ICA.BAD.length = 4;
    cfg.STIM.ICA.BAD.percentbad = 70;
    cfg.STIM.ICA.BAD.freqlim50 = 0;
    cfg.STIM.ICA.BAD.freqlim60 = 0;
    cfg.STIM.ICA.BAD.freqlimlow = 0;
    cfg.STIM.ICA.BAD.freqlimhigh = 0;
    cfg.STIM.ICA.BAD.spectslow = [0.5 2];
    cfg.STIM.ICA.BAD.spectshigh = [15 50];
    cfg.STIM.ICA.BAD.zvaluebroken = 0;
    cfg.STIM.ICA.BAD.zvaluevars = 0;
    cfg.STIM.ICA.BAD.zvaluehurst = 0;
    cfg.STIM.ICA.BAD.zvaluemedian = 0;
    cfg.STIM.ICA.BAD.zvalueamplitude = 0;
    cfg.STIM.ICA.BAD.zvaluekurtosis = 0;
    cfg.STIM.ICA.BAD.zvalueeye = 0;
    cfg.STIM.ICA.BAD.eog = [];
    cfg.STIM.ICA.BAD.ecgdetect = 3;
    cfg.STIM.ICA.BAD.PEAK2MIN = [];
    cfg.STIM.ICA.BAD.AVG = [];
    cfg.STIM.ICA.BAD.AVG.marker = {'STIM'};
    cfg.STIM.ICA.BAD.AVG.combinemarker = true;
    if isfield(header,'samplingrate')
        cfg.STIM.ICA.BAD.AVG.markerlength = round(header.samplingrate/130) + 1;
    else
        cfg.STIM.ICA.BAD.AVG.markerlength = 9;
    end
    cfg.STIM.ICA.BAD.AVG.markerOffset = floor(cfg.STIM.ICA.BAD.AVG.markerlength/2);
    cfg.STIM.ICA.BAD.AVG.AVGmethod = 'mean';
    cfg.STIM.ICA.BAD.AVG.Reject = [];
    cfg.STIM.ICA.BAD.AVG.correctbaseline = 'disabled';
    cfg.STIM.ICA.BAD.AVGmode = 'detect bad';
    cfg.STIM.ICA.BAD.AVGstd = 2.5;
    cfg.STIM.ICA.ICABACK = [];
end

Prompt = cell(0,2);
Formats = [];

Prompt(end+1,:) = {'Output-folder', 'folder'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'text';
Formats(end,1).size = 100;

Prompt(end+1,:) = {'Detect DBS stimulator','DETECT'};
Formats(end+1,1).type = 'check';
Formats(end,1).callback = {@detect_dbs,'DETECT','DETECT'};

Prompt(end+1,:) = {'Remove principal component','PCA'};
Formats(end+1,1).type = 'check';

Prompt(end+1,:) = {'Independent component analysis','ICA'};
Formats(end+1,1).type = 'check';
Formats(end,1).callback = {@do_ICA,'ICA','ICA',header};

[cfg.STIM,Cancelled] = inputsdlg(Prompt,'Remove DBS',Formats,cfg.STIM);
pause(0.1);
if isempty(cfg.STIM) | Cancelled == 1
    cfg.STIM = [];
    skipprocessing = 1;
end

end

function settings = detect_dbs(settings)
    if ~isfield(settings,'stimrange') | isempty(settings.stimrange)
        settings.stimrange = [100,150;80,100;150,250];
    end
    if ~isfield(settings,'stimmaxpercent') | isempty(settings.stimmaxpercent)
        settings.stimmaxpercent = 70;
    end
    
    Prompt = cell(0,2);
    Formats = [];
    
    Prompt(end+1,:) = {'   Range (Hz)','stimrange'};
    Formats(end+1,1).type = 'edit';
    Formats(end,1).format = 'result';
    Formats(end,1).size = [70 20];
    Formats(end,1).callback = {@lab_table_dialog,'stimrange', ...
        'stimrange',{'lowfreq','highfreq'},'Range (Hz)',1};
    
    Prompt(end+1,:) = {'   Min amplitude for first range (percent)','stimmaxpercent'};
    Formats(end+1,1).type = 'edit';
    Formats(end,1).format = 'float';
    Formats(end,1).size = [30];
    Formats(end,1).limits = [0 100];
    
    [settings,Cancelled] = inputsdlg(Prompt,'Detect DBS',Formats,settings);
    pause(0.1);
    if isempty(settings) | Cancelled == 1
        settings = [];
    end
end

function settings = do_ICA(settings,header)
    if isempty(settings)
        settings.type = 'runica';
        settings.extended = 0;
        settings.pca = 0;
        settings.good = true;
        settings.epoch = 0;
        settings.plottopo = false;
        settings.BAD = [];
        settings.BAD.length = 4;
        settings.BAD.percentbad = 70;
        settings.BAD.freqlim50 = 0;
        settings.BAD.freqlim60 = 0;
        settings.BAD.freqlimlow = 0;
        settings.BAD.freqlimhigh = 0;
        settings.BAD.spectslow = [0.5 2];
        settings.BAD.spectshigh = [15 50];
        settings.BAD.zvaluebroken = 0;
        settings.BAD.zvaluevars = 0;
        settings.BAD.zvaluehurst = 0;
        settings.BAD.zvaluemedian = 0;
        settings.BAD.zvalueamplitude = 0;
        settings.BAD.zvaluekurtosis = 0;
        settings.BAD.zvalueeye = 0;
        settings.BAD.eog = [];
        settings.BAD.ecgdetect = 3;
        settings.BAD.PEAK2MIN = [];
        settings.BAD.AVG = [];
        settings.BAD.AVG.marker = {'STIM'};
        settings.BAD.AVG.combinemarker = true;
        if isfield(header,'samplingrate')
            settings.BAD.AVG.markerlength = round(header.samplingrate/130) + 1;
        else
            settings.BAD.AVG.markerlength = 9;
        end
        settings.BAD.AVG.markerOffset = floor(settings.BAD.AVG.markerlength/2);
        settings.BAD.AVG.AVGmethod = 'mean';
        settings.BAD.AVG.Reject = [];
        settings.BAD.AVG.correctbaseline = 'disabled';
        settings.BAD.AVGmode = 'detect bad';
        settings.BAD.AVGstd = 2.5;
        settings.ICABACK = [];
    end
    Prompt = cell(0,2);
    Formats = [];
        
    Prompt(end+1,:) = {'Std of AVG >','AVGstd'};
    Formats(end+1,1).type = 'edit';
    Formats(end,1).format = 'float';
    Formats(end,1).size = 30;
    Formats(end,1).limits = [0 100];
    
    [settings.BAD,Cancelled] = inputsdlg(Prompt,'ICA - Detect DBS',Formats,settings.BAD);
    pause(0.1);
    if isempty(settings.BAD) | Cancelled == 1
        settings = [];
    end
end