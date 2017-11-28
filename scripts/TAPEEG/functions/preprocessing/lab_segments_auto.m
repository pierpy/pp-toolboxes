% Script for automated segmentation of eeg/meg, called by
% lab_define_segments
%
% [Segments,events] = lab_segments_auto(data,header,cfg)
%
% written by F. Hatz 2013

function [Segments,events] = lab_segments_auto(data,header,cfg)
disp('   Automated segments selection')

if ~isfield(cfg,'SEG') | ~isfield(cfg.SEG,'AUTO') | isempty(cfg.SEG.AUTO)
    if ~isfield(cfg,'MAIN') | ~isfield(cfg.MAIN,'auto') | cfg.MAIN.auto == 0
        [cfg.SEG.AUTO,skipprocessing] = lab_set_segments_auto([],header);
        if skipprocessing == 1
            Segments = [1 header.numtimeframes];
            return
        end
    else
        Segments = [1 header.numtimeframes];
        return
    end
end

% filtering with butterworth filter
if isfield(cfg,'FILT') & ~isempty(cfg.FILT) & isfield(cfg.FILT,'Filter') & ~isempty(cfg.FILT.Filter)
    cfg.SEG.AUTO.FILT = cfg.FILT;
    cfg.SEG.AUTO.FILT.filtermode = 'butter';
    cfg.SEG.AUTO.FILT.filtorder = 2;
else
    cfg.SEG.AUTO.FILT.filtermode = 'butter';
    cfg.SEG.AUTO.FILT.filtorder = 2;
    cfg.SEG.AUTO.FILT.Filter.highpass = 1;
    cfg.SEG.AUTO.FILT.Filter.lowpass = 70;
    cfg.SEG.AUTO.FILT.Filter.notch = [];
    cfg.SEG.AUTO.FILT.filterauxchannels = false;
end
cfg.SEG.AUTO.FILT.nodisp = 1;
[data,header,cfg.SEG.AUTO] = lab_filtering(data,header,cfg.SEG.AUTO,1);

% correct markers
cfg.SEG.AUTO = lab_markerselect(cfg.SEG.AUTO,cfg,header);

% Interpolate bad channels for analysis
if ~isfield(header,'goodchans')
    header.goodchans = 1:header.numdatachannels;
    header.badchans = [];
end

% Find epochs
if size(data,2) < cfg.SEG.AUTO.minduration*header.samplingrate
    Segments = [1 header.numtimeframes];
    return
end

settings = cfg.SEG.AUTO.BAD;
stepsize = settings.length;
shiftnr = ceil(cfg.SEG.AUTO.minduration / stepsize);
shift = header.samplingrate * stepsize;

[~,bad] = lab_detect_bad(data,header,settings,[],'noverbose');
electrodesValid = abs(bad.epochs - 1);

electrodesMarker = ones(size(electrodesValid));
% search bad shifts by markers
if isfield(settings,'markerexclude') & ~isempty(settings.markerexclude) & isfield(header,'events')
    startE = [];
    stopE = [];
    for i = 1:size(settings.markerexclude,1)
        tmp = find(strcmp(header.events.TYP,settings.markerexclude(i,1)));
        if ~isempty(tmp)
            startE = [startE header.events.POS(tmp)]; %#ok<AGROW>
            stopE = [stopE (header.events.POS(tmp)+header.events.DUR(tmp))]; %#ok<AGROW>
        end
    end
    startE = ceil(double(startE) / shift);
    stopE = ceil(double(stopE) / shift);
    for i = 1:size(startE,2)
        electrodesValid(:,startE(i):stopE(i)) = 0;
        electrodesMarker(:,startE(i):stopE(i)) = 0;
    end
    clearvars startE stopE
end

% Look for marker start and label all epochs before that marker as bad
if isfield(settings,'markerstart') & ~isempty(settings.markerstart) & isfield(header,'events')
    tmp = find(strcmp(header.events.TYP,settings.markerstart), 1, 'last' );
    tmp = header.events.POS(tmp);
    if ~isempty(tmp)
        electrodesValid(:,1:ceil(tmp / shift)) = 0;
    end
end

% Look for marker stop and label all epochs after that marker as bad
if isfield(settings,'markerstop') & ~isempty(settings.markerstop) & isfield(header,'events')
    tmp = find(strcmp(header.events.TYP,settings.markerstop), 1 );
    tmp = header.events.POS(tmp);
    if ~isempty(tmp)
        electrodesValid(:,ceil(tmp/shift):end) = 0;
    end
end

% restrict good epochs to markers
if isfield(settings,'markerinclude') & ~isempty(settings.markerinclude) & isfield(header,'events')
    electrodesValid2 = zeros(size(electrodesValid));
    startE = [];
    stopE = [];
    for i = 1:size(settings.markerinclude,1)
        tmp = find(strcmp(header.events.TYP,settings.markerinclude(i,1)));
        if ~isempty(tmp)
            startE = [startE header.events.POS(tmp)]; %#ok<AGROW>
            stopE = [stopE (header.events.POS(tmp)+header.events.DUR(tmp))]; %#ok<AGROW>
        end
        startE = ceil(double(startE) / shift);
        stopE = floor(double(stopE) / shift);
    end
    for i = 1:size(startE,2)
        electrodesValid2(:,startE(i):stopE(i)) = 1;
    end
    electrodesValid(~electrodesValid2) = 0;
    clearvars electrodesValid2 startE stopE
end

% set minimalpart to integer value = number of epochs
minimalpart = ceil(cfg.SEG.AUTO.durationall / stepsize);
if minimalpart*shift > size(data,2)
    disp('    Duration of all segments to long, taking all input file')
    Segments = [1 header.numtimeframes];
    return
end

% evaluate peak2min
electrodesValid = mean(electrodesValid,1);
if isfield(bad,'peak2min') & ~isempty(bad.peak2min)
    electrodesValid = electrodesValid .* bad.peak2min;
end

% find good epochs (sensitivity is increased until minimalpart is reached)
if isfield(cfg.SEG.AUTO,'percentgood')
    percentgood = cfg.SEG.AUTO.percentgood / 100;
else
    percentgood = 1;
end
goodSeg = [];
while sum(goodSeg) < minimalpart
    goodSeg = (electrodesValid >= percentgood);
    badSeg = find(electrodesValid < percentgood);
    tmp = 0;
    for i = 1:size(badSeg,2)
        if (badSeg(i) - tmp) < shiftnr
            goodSeg((tmp+1):badSeg(i)) = false;
        end
        if (badSeg(i) + shiftnr - 1) > size(goodSeg,2)
            goodSeg(badSeg(i):end) = false;
        end
        tmp = badSeg(i);
    end
    if sum(goodSeg) < minimalpart
        percentgood = percentgood * 0.99;
    end
end

SegStart = find(diff(goodSeg,1)==1)+1;
SegStop = find(diff(goodSeg,1)==-1);
if goodSeg(1) == 1;
    SegStart = [1 SegStart];
end
if goodSeg(end) == 1
    SegStop = [SegStop size(goodSeg,2)];
end

% redcue number of segments if summarized length > durationall
allpart = 0;
nseg = 0;
while allpart<minimalpart
    nseg = nseg+1;
    allpart = sum((SegStop(1:nseg) - SegStart(1:nseg))+1);
end
SegStart = SegStart(1:nseg);
SegStop = SegStop(1:nseg);
clearvars allpart nseg

if length(SegStart) ~= length(SegStop)
    Segments = [1 header.numtimeframes];
    return
else
    SegStart = (SegStart-1) * shift + 1;
    SegStop = SegStop * shift;
    Segments = [SegStart' SegStop'];
end

for i = 1:size(Segments,1)
    events.POS(1,(i-1)*2+1) = int64(Segments(i,1));
    events.DUR(1,(i-1)*2+1) = int64(1);
    events.OFF(1,(i-1)*2+1) = int64(0);
    events.TYP{1,(i-1)*2+1} = ['AutoSegStart' num2str(i)];
    events.POS(1,(i-1)*2+2) = int64(Segments(i,2));
    events.DUR(1,(i-1)*2+2) = int64(1);
    events.OFF(1,(i-1)*2+2) = int64(0);
    events.TYP{1,(i-1)*2+2} = ['AutoSegStop' num2str(i)];
end

% write edf if necessary
if isfield(cfg.SEG.AUTO,'montage') & ~isempty(cfg.SEG.AUTO.montage) & isfield(cfg,'EEG_file') & isfield(cfg,'EEG_filepath')
    header = lab_mix_markers(header,events);
    [~,~,~,EDF_file] = lab_filename(cfg.EEG_file);
    disp('   Write *AutoSeg.edf')
    cfgtmp.EDF.eegsource = 'montage';
    cfgtmp.EDF.interpolatebad = false;
    cfgtmp.EDF.filepath = '';
    cfgtmp.EDF.mountpath = '';
    cfgtmp.EDF.montage = cfg.SEG.AUTO.montage;
    cfgtmp.EDF_file = [EDF_file '_AutoSeg.edf'];
    cfgtmp.EDF_filepath = cfg.EEG_filepath;
    [~,~,cfgtmp] = lab_export2edf(data,header,cfgtmp,'nofolder');
    cfg.SEG.AUTO.montage = cfgtmp.EDF.montage;
end
%--------------------------------------------------------------------------
% Write verbose file (*.vrb)
%--------------------------------------------------------------------------
if isfield(cfg,'EEG_file') & isfield(cfg,'EEG_filepath')
    [~,~,~,Verbose_file] = lab_filename(cfg.EEG_file);
    Verbose_file = fullfile(cfg.EEG_filepath,[Verbose_file '_AutoSeg.vrb']);
    fid=fopen( Verbose_file,'w');
    fprintf(fid,'Auto Segments\n');
    fprintf(fid,['EEG-file: ' cfg.EEG_file]);
    fprintf(fid,'\n\n');
    fprintf(fid,['Minimal length of segments (seconds): ' num2str(cfg.SEG.AUTO.minduration)]);
    fprintf(fid,'\n');
    fprintf(fid,['Summed length of segments (seconds): ' num2str(sum(goodSeg)*2)]);
    fprintf(fid,'\n');
    if isfield(cfg.SEG.AUTO,'BAD')
        fprintf(fid,'Settings for detecting bad channels:\n');
        lab_write_bad_vrb(fid,cfg.SEG.AUTO.BAD,bad);
        fprintf(fid,'\n\n');
    end
    fprintf(fid,['Percent of minimal good channels in epoch: ' num2str(percentgood*100)]);
    fprintf(fid,'\n');
    if ~isempty(cfg.SEG.AUTO.markerexclude)
        fprintf(fid,['Markers excluded: ' sprintf('%s|',cfg.SEG.AUTO.markerexclude{:})]);
        fprintf(fid,'\n');
    end
    if ~isempty(cfg.SEG.AUTO.markerinclude)
        fprintf(fid,['Markers included: ' sprintf('%s|',cfg.SEG.AUTO.markerinclude{:})]);
        fprintf(fid,'\n');
    end
    fclose(fid);
    
    if ~isempty(cfg.SEG.AUTO.markerexclude) & strcmp(cfg.SEG.AUTO.markerexclude{1,1},'all')
        cfg.SEG.AUTO.markerexclude = cellstr('all');
    end
    if ~isempty(cfg.SEG.AUTO.markerinclude) & strcmp(cfg.SEG.AUTO.markerinclude{1,1},'all')
        cfg.SEG.AUTO.markerinclude = cellstr('all');
    end
end

return