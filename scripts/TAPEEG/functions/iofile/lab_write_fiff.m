% Write Elekty .fif (MNE toolbox)
%
% lab_write_fiff(filename,data,header)
%
% written by F. Hatz 2012

function lab_write_fiff(filename,data,header,cfg)
[~,filepath,~,filename] = lab_filename(filename);
filename = fullfile(filepath,[filename '.fif']);

skipprocessing = 0;
if ~exist('cfg','var')
    cfg = [];
end

if isfield(header,'orig') & isfield(header.orig,'info')
    header.orig.sfreq = header.samplingrate;
    header.orig.last_samp = header.orig.first_samp + header.numtimeframes -1;
    if isfield(header,'includechans')
        header.orig.info.chs = header.orig.info.chs(1,header.includechans);
        header.orig.info.ch_names = header.orig.info.ch_names(1,header.includechans);
        header.orig.info.nchan = header.numchannels;
    end
    if isfield(header,'lowpass')
        header.orig.info.lowpass = header.lowpass;
        header.orig.info.highpass = header.highpass;
    end
else
    disp('   Abort - no FIF-information in header')
    skipprocessing = 1;
end

if isfield(cfg,'SCALE') & isfield(cfg.SCALE,'scales')
    for i = 1:size(cfg.SCALE.scales,1)
        data(cfg.SCALE.scales(i,1):cfg.SCALE.scales(i,2),:) = data(cfg.SCALE.scales(i,1):cfg.SCALE.scales(i,2),:) / cfg.SCALE.scales(i,3);
    end
end

if skipprocessing == 0
    FIFF = fiff_define_constants();
    
    % select only datachannels (other channels not supported yet)
    want_meg   = true;
    want_eeg   = false;
    want_stim  = false;
    picks = fiff_pick_types(header.orig.info,want_meg,want_eeg,want_stim);
    
    % Match channels to fiff-header information
    channelstmp = header.orig.info.ch_names;
    channelstmp = channelstmp(1,picks);
    channelstmp2 = cellstr(header.channels)';
    resort = zeros(1,size(channelstmp,2));
    for i = 1:size(channelstmp,2)
        resort(1,i) = find(strcmp(channelstmp2,channelstmp(1,i))==1);
    end
    if  find(resort == 0) > 0
        tmp = find(resort > 0);
        if ~isempty(tmp)
            disp('   Warning: some original datachannels missing')
            picks = picks(tmp);
            resort = resort(tmp);
        else
            disp('   Abort: datachannels missing')
            skipprocessing = 1;
        end
        clearvars tmp
    end
    clearvars channelstmp channelstmp2
end

if skipprocessing == 0
    % resort data to match fiff-header
    data = data(resort,:);
    clearvars resort
    
    % Write data
    [outfid,cals] = fiff_start_writing_raw(filename,header.orig.info,picks);
    
    from        = header.orig.first_samp;
    to          = header.orig.last_samp;
    %quantum_sec = 10;
    %quantum     = ceil(quantum_sec*header.orig.info.sfreq);
    quantum     = to - from + 1; % write data in one block
    
    first_buffer = true;
    for first = from:quantum:to
        last = first+quantum-1;
        if last > to
            last = to;
        end
        if first_buffer
            if first > 0
                fiff_write_int(outfid,FIFF.FIFF_FIRST_SAMPLE,first);
            end
            first_buffer = false;
        end
        fiff_write_raw_buffer(outfid,data,cals);
    end
    fiff_finish_writing_raw(outfid);
    
    % Write EEGinfo-file (*.txt)
    lab_write_eeginfo(filename,header)
end