% Resample EEG
%
% [data,header,cfg]=lab_resample_data(data,header,cfg)
%
%   data                = data (chans x timeframes)
%   header.samplingrate = actual sampling rate
%   new_samplingrate    = new EEG rate
%
% Author:  Florian Hatz, Neurology Basel 2012

function [data,header,cfg,skipprocessing]=lab_resample_data(data,header,cfg,verbosefile)

skipprocessing = 0;

if exist('cfg','var') & isnumeric(cfg)
    tmp =cfg;
    cfg = [];
    cfg.RESAMPLE.new_samplingrate = tmp;
    clearvars tmp;
end

if exist('header','var') & ~isempty(header) & isstruct(header) & isfield(header,'samplingrate')
    old_samplingrate = header.samplingrate;
elseif exist('header','var') & ~isempty(header) & isnumeric(header)
    old_samplingrate = header;
    header = [];
else
    disp ('   Ask for old and new samplingrate')
    prompt={'old samplingrate','new samplingrate'};
    name='Settings resampling';
    numlines = [1,40;1,40];
    answer=inputdlg(prompt,name,numlines);
    pause(0.1)
    if size(answer,1) > 0
        old_samplingrate = str2num(answer{1,1}); %#ok<ST2NM>
        cfg.RESAMPLE.new_samplingrate = str2num(answer{2,1}); %#ok<ST2NM>
    else
        skipprocessing = 1;
        old_samplingrate = [];
        cfg.RESAMPLE.new_samplingrate=[];
    end
end

if ~exist('cfg','var') | ~isfield(cfg,'RESAMPLE') | ~isfield(cfg.RESAMPLE,'new_samplingrate') | ...
        isempty(cfg.RESAMPLE.new_samplingrate)
    disp ('   Ask for new samplingrate')
    prompt={'new samplingrate'};
    name='Settings resampling';
    numlines(1,1) = 1;
    numlines(1,2) = 40;
    defaultanswer={num2str(old_samplingrate)};
    answer=inputdlg(prompt,name,numlines,defaultanswer);
    pause(0.1)
    if size(answer,1) > 0
        cfg.RESAMPLE.new_samplingrate = str2num(answer{1,1}); %#ok<ST2NM>
    else
        skipprocessing = 1;
        cfg.RESAMPLE.new_samplingrate=[];
    end
end

if ~isempty(cfg.RESAMPLE.new_samplingrate) & ~isempty(old_samplingrate) & ...
        old_samplingrate ~= cfg.RESAMPLE.new_samplingrate
    if exist('verbosefile','var')
        disp(['   Resampling to ' num2str(cfg.RESAMPLE.new_samplingrate) 'Hz'])
    end
    numchannels = size(data,1);
    lengthin = size(data,2);
    [N,D] = rat(cfg.RESAMPLE.new_samplingrate/old_samplingrate,0.0001);
    lengthout = floor(lengthin*N/D);
    lengthin = floor(lengthout*D/N);
    clearvars tmp
    for trial = 1:size(data,3);
        for ch = 1:numchannels
            tmp(ch,:,trial) = resample(data(ch,1:lengthin,trial),N,D);
        end
    end
    data = tmp;
    clearvars tmp
    
    header.samplingrate = cfg.RESAMPLE.new_samplingrate;
    header.numtimeframes = size(data,2);
    if isfield(header,'events')
        header.events.POS = int64(round(header.events.POS ./ (D/N)));
        header.events.POS(header.events.POS == 0) = 1;
        header.events.DUR = int64(ceil(header.events.DUR ./ (D/N)));
        header.events.DUR(header.events.DUR < 1) = int64(1);
        if isfield(header.events,'OFF')
            header.events.OFF = int64(round(header.events.OFF ./ (D/N)));
        end
    end
    
    if exist('verbosefile','var')
        fid=fopen(verbosefile,'w');
        fprintf(fid,'Re-sampling data\n');
        fprintf(fid,datestr(now,0));
        fprintf(fid,'\n');
        fprintf(fid,['Input file: ' header.EEG_file]);
        fprintf(fid,'\n\n');
        fprintf(fid,['New samplingrate: ' num2str(cfg.RESAMPLE.new_samplingrate) 'Hz']);
        fprintf(fid,'\n');
        fclose(fid);
    end
else
    disp('   Samplingrate unchanged, no resampling needed')
end

return;