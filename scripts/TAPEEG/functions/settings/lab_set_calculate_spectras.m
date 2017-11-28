function [cfg,skipprocessing] = lab_set_calculate_spectras(cfg,header,doshort)

skipprocessing = 0;

if ~exist('doshort','var')
    doshort = 0;
end

global THeader
if ~isempty(THeader)
    header = THeader;
elseif ~exist('header','var')
    header = [];
end

if isfield(header,'numdatachannels')
    numchans = header.numdatachannels;
else
    numchans = [];
end
if doshort == 1
    tmp = cfg;
    cfg = [];
    cfg.FFT = tmp;
    clearvars tmp
end
if ~exist('cfg','var') | isempty(cfg) | ~isfield(cfg,'FFT') | ~isfield(cfg.FFT,'winsize')
    if doshort == 0
        cfg.FFT.folder = 'Spect';
        cfg.FFT.eegsource = 'input';
        cfg.FFT.method = 'multitaper';
        cfg.FFT.deleteold = true;
        cfg.FFT.interpolatemethod = 'spherical';
        cfg.FFT.percentgood = 90;
        cfg.FFT.Mappings = [];
        cfg.FFT.plotspectras = false;
        cfg.FFT.hanpercent = 0;
        cfg.FFT.calcpower = false;
    else
        cfg.FFT.hanpercent = 20;
    end
    cfg.FFT.winsize = 4;
    cfg.FFT.lowfreq = 1;
    cfg.FFT.highfreq = 50;
    cfg.FFT.spectralbands = lab_get_spectralbands;
    cfg.FFT.spectralbandsI = false;
    cfg.FFT.BAD.freqlim50 = 50;
    cfg.FFT.BAD.freqlim60 = 50;
    cfg.FFT.BAD.freqlimlow = 70;
    cfg.FFT.BAD.freqlimhigh = 50;
    cfg.FFT.BAD.spectslow = [0.5 2];
    cfg.FFT.BAD.spectshigh = [15 50];
    cfg.FFT.BAD.zvaluebroken = 0;
    if isfield(header,'interpolated') & ~isempty(header.interpolated)
        cfg.FFT.BAD.zvaluevars = 4;
    else
        cfg.FFT.BAD.zvaluevars = 3;
    end
    cfg.FFT.BAD.zvaluehurst = 3;
    cfg.FFT.BAD.zvaluemedian = 3;
    cfg.FFT.BAD.zvaluecorr = 0;
    cfg.FFT.BAD.LAPL = [];
    cfg.FFT.BAD.markerexclude = {};
end
if isfield(cfg.FFT,'spectralbands') & isnumeric(cfg.FFT.spectralbands) & size(cfg.FFT.spectralbands,2) == 2
    cfg.FFT.spectralbands = cat(2,cell(size(cfg.FFT.spectralbands,1),1),num2cell(cfg.FFT.spectralbands));
end

Prompt = cell(0,2);
Formats = [];
if doshort == 0
    Prompt(end+1,:) = {'Output-folder', 'folder'};
    Formats(end+1,1).type = 'edit';
    Formats(end,1).format = 'text';
    Formats(end,1).size = 100;
    Formats(end,1).span = [1 3];
    
    Prompt(end+1,:) = {'Reference','eegsource'};
    Formats(end+1,1).type = 'list';
    Formats(end,1).style = 'popupmenu';
    Formats(end,1).format = 'input';
    if ~max(strcmp({'channels','mean','median','laplacian','montage','input'},cfg.FFT.eegsource))
        Formats(end,1).items = {cfg.FFT.eegsource,'channels','mean','median','laplacian','montage','input'};
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

    Prompt(end+1,:) = {'Interpolate bad channels (method)','interpolatemethod'};
    Formats(end+1,1).type = 'list';
    Formats(end,1).style = 'popupmenu';
    Formats(end,1).format = 'input';
    Formats(end,1).items = {'disabled';'spherical';'3D'};
    Formats(end,1).span = [1 3];
    
    Formats(end+1,1).type = 'none';
    Formats(end,1).span = [1 3];
    
    Prompt(end+1,:) = {'Detect bad periods','BAD'};
    Formats(end+1,1).type = 'edit';
    Formats(end,1).format = 'result';
    Formats(end,1).callback = {@lab_set_detect_bad,'BAD','BAD',[],[],0,0,0,0,0,0,0,0,1,0};
    Formats(end,1).size = [150 185];
    Formats(end,1).span = [1 3];
    
    Prompt(end+1,:) = {'Minimum percent good channels in single window','percentgood'};
    Formats(end+1,1).type = 'edit';
    Formats(end,1).format = 'float';
    Formats(end,1).limits = [0 100];
    Formats(end,1).size = 35;

    Formats(end+1,1).type = 'none';
    Formats(end,1).span = [1 3];
    
    Prompt(end+1,:) = {'Method','method'};
    Formats(end+1,1).type = 'list';
    Formats(end,1).style = 'popupmenu';
    Formats(end,1).format = 'input';
    Formats(end,1).items = {'welch','multitaper','wavelet','rihaczek'};
    Formats(end,1).span = [1 3];
    Formats(end,1).callback = {@set_method,'@ALL','@ALL'};
end

Prompt(end+1,:) = {'FFT window (seconds)','winsize'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'float';
Formats(end,1).limits = [0 inf];
Formats(end,1).size = 35;

Prompt(end+1,:) = {'Hanning-window (percent)','hanpercent'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'float';
Formats(end,1).limits = [0 100];
Formats(end,1).size = 35;
Formats(end,1).span = [1 2];

Formats(end+1,1).type = 'none';
Formats(end,1).span = [1 3];

Prompt(end+1,:) = {'Lowest frequency','lowfreq'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'float';
Formats(end,1).limits = [0 inf];
Formats(end,1).size = 35;

Prompt(end+1,:) = {'Highest frequency','highfreq'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'float';
Formats(end,1).limits = [0 inf];
Formats(end,1).size = 35;
Formats(end,1).span = [1 2];

Prompt(end+1,:) = {'Spectral Bands','spectralbands'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'result';
Formats(end,1).callback = {@lab_get_spectralbands,'spectralbands','spectralbands','spectralbandsI',true};
Formats(end,1).span = [1 2];

Prompt(end+1,:) = {'Individual Bands','spectralbandsI'};
Formats(end+1,1).type = 'check';
Formats(end,1).format = 'input';
Formats(end,1).callback = {@do_spectralbandsI,'@ALL','@ALL'};

if doshort == 0
    Prompt(end+1,:) = {'Calculate power per channel','calcpower'};
    Formats(end+1,1).type = 'check';
    
    Prompt(end+1,:) = {'Reduce results to Mappings','Mappings'};
    Formats(end+1,1).type = 'check';
    Formats(end,1).callback = {@lab_load_mappings,'Mappings','Mappings',cfg,numchans};
    Formats(end,1).span = [1 2];
    
    Formats(end+1,1).type = 'none';
    Formats(end,1).span = [1 3];
    
    Prompt(end+1,:) = {'Write Excel-files','writexls'};
    Formats(end+1,1).type = 'check';
    
    Prompt(end+1,:) = {'Plot spectras','plotspectras'};
    Formats(end+1,1).type = 'check';
    Formats(end,1).span = [1 2];
    
    Prompt(end+1,:) = {'Delete results from previous run','deleteold'};
    Formats(end+1,1).type = 'check';
    Formats(end,1).span = [1 3];
end

[cfg.FFT,Cancelled] = inputsdlg(Prompt,'Spectral Analysis',Formats,cfg.FFT);
if isempty(cfg.FFT) | Cancelled == 1
    cfg.FFT = [];
    skipprocessing = 1;
    return
end

if doshort == 1
    cfg = cfg.FFT;
end

end

function settings = set_method(settings)
   if strcmp(settings.method,'multitaper')
       settings.hanpercent = 0;
   else
       settings.hanpercent = 0.2;
   end
end

function settings = do_spectralbandsI(settings)
    if settings.spectralbandsI == false
        settings.spectralbandsI = true;
        settings.spectralbands = [1,2,4,5,6,8,9];
    else
        settings.spectralbandsI = false;
        settings.spectralbands = lab_get_spectralbands;
    end
end