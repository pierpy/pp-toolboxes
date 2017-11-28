function [data,header,cfg] = lab_detect_eyeblinks(data,header,cfg)

if ~exist('header','var') | ~isfield(header,'samplingrate')
    header = lab_create_header(data);
end
if ~exist('cfg','var')
    cfg = [];
end

if ~isfield(cfg,'EEG_file')
    cfg.EEG_file = header.EEG_file;
    cfg.EEG_filepath = header.EEG_filepath;
end
[~,~,~,cfg.EEG_fileS] = lab_filename(cfg.EEG_file);

if ~isfield(header,'goodchans') | isempty(header.goodchans)
    header.goodchans = 1:header.numdatachannels;
    header.badchans = [];
end

if ~isfield(cfg,'EYEBLINKS') | ~isfield(cfg.EYEBLINKS,'eog')
    [cfg,skipprocessing] = lab_set_detect_eyeblinks(cfg,header);
    if skipprocessing == 1
        return
    end
end

[datatmp,headertmp] = lab_interpolate_bad(data,header);
settings = lab_calculate_eog(datatmp,headertmp,cfg.EYEBLINKS);
blinkchans = settings.blinkchans;
clearvars datatmp headertmp settings
Eyeblinks = zeros(1,size(data,2));
for i = 1:size(blinkchans,1)
    tmp = find(blinkchans(i,:) > cfg.EYEBLINKS.SD * std(blinkchans(i,:)));
    if ~isempty(tmp)
        Eyeblinks(tmp) = Eyeblinks(tmp) + 1;
    end
end

if ~isfield(cfg.EYEBLINKS,'MinChans') | isempty(cfg.EYEBLINKS.MinChans)
    cfg.EYEBLINKS.MinChans = 66;
end
Eyeblinks = find(Eyeblinks >= (cfg.EYEBLINKS.MinChans*size(blinkchans,1)/100));
if ~isempty(Eyeblinks)
    tmp = find(diff(Eyeblinks) > 1);
    tmp = [0 tmp length(Eyeblinks)];
    tmp2 = [];
    for i = 1:length(tmp)-1
        tmp2(end+1) = round(mean(Eyeblinks(tmp(i)+1:tmp(i+1)))); %#ok<AGROW>
    end
    Eyeblinks = tmp2;
    clearvars tmp tmp2
    
    Eflag = zeros(1,size(data,2));
    for i = 1:length(Eyeblinks)
        StartE = Eyeblinks(i) - round((cfg.EYEBLINKS.duration/2) * header.samplingrate);
        StopE = Eyeblinks(i) + round((cfg.EYEBLINKS.duration/2) * header.samplingrate);
        if StartE < 1
            StartE = 1;
        end
        if StopE > size(data,2)
            StopE = size(data,2);
        end
        Eflag(StartE:StopE) = 1;
    end
    
    events = [];
    tmp = find(abs(diff(Eflag))>0);
    tmp = tmp(:)';
    EndEflag = [tmp size(data,2)];
    StartEflag = [1 tmp+1];
    Idx = find(Eflag(StartEflag) ~= 0);
    if ~isempty(Idx)
        StartEflag = StartEflag(Idx);
        EndEflag = EndEflag(Idx);
        events.POS = int64(StartEflag(:)');
        events.DUR = int64(EndEflag(:)' - StartEflag(:)' + 1);
        events.OFF = int64(zeros(1,length(StartEflag)));
        events.TYP = repmat(cellstr('eyeblink'),1,length(StartEflag));
    end
    if ~isempty(events)
        if isfield(header,'events') & ~isempty(header.events)
            header = lab_mix_markers(header,events);
        else
            header.events = events;
        end
    end
else
    disp('     no eyeblinks detected')
end

end