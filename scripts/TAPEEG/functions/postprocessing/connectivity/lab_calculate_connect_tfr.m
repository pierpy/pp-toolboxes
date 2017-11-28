% Calculate connectivity based on Rihaczek time-frequency distribution in
% specified frequency ranges
%
% [result,cfg] = lab_calculate_connect_tfr(data,header,cfg)
%
% written by F. Hatz 2012

function [Result,cfg,header] = lab_calculate_connect_tfr(data,header,cfg)

if isfield(cfg,'Output_file') & ~isempty(cfg.Output_file)
    [~,~,~,cfg.Output_fileS] = lab_filename(cfg.Output_file);
end

if ~isfield(cfg.CONNECT,'PHASE') | ~isfield(cfg.CONNECT.PHASE,'window')
    cfg.CONNECT.PHASE.window = 4;
end
if ~isfield(cfg.CONNECT.PHASE,'freqwindow')
    cfg.CONNECT.PHASE.freqwindow = cfg.CONNECT.PHASE.window * header.samplingrate;
end

if ~isfield(cfg.CONNECT,'freqs')
    cfg.CONNECT.freqs = [1 4;4 8;8 10;10 13;13 30];
end
if size(cfg.CONNECT.freqs,2) == 4
    Freqs = cfg.CONNECT.freqs;
elseif size(cfg.CONNECT.freqs,2) == 2
    Freqs = [cfg.CONNECT.freqs cfg.CONNECT.freqs];
else
    Freqs = [0 0 0 0];
end

if ~isfield(cfg.CONNECT.PHASE,'step')
    cfg.CONNECT.PHASE.step = 8;
end

% set temporary variables
Dindex = 1:size(data,2);
datatmp = data;
headertmp = header;

% prepare exclude markers if necessary
if isfield(cfg.CONNECT,'MARKER') & isfield(cfg.CONNECT.MARKER,'markerexclude') & ...
        iscell(cfg.CONNECT.MARKER.markerexclude) & ~isempty(cfg.CONNECT.MARKER.markerexclude) & ...
        isfield(headertmp,'events') & ~isempty(headertmp.events)
    disp('    Exclude periods with selected markers');
    [datatmp,headertmp,Idx2] = lab_exclude_markers(datatmp,headertmp,cfg.CONNECT.MARKER.markerexclude);
    if isempty(datatmp)
        disp('     Abort: calculation of connectivity not possible');
        return
    end
    Dindex = Dindex(1,Idx2);
    clearvars Idx2
end

% prepare include markers if necessary
if isfield(cfg.CONNECT,'MARKER') & isfield(cfg.CONNECT.MARKER,'markerinclude') & ...
        iscell(cfg.CONNECT.MARKER.markerinclude) & ~isempty(cfg.CONNECT.MARKER.markerinclude) & ...
        isfield(headertmp,'events') & ~isempty(headertmp.events)
    disp('    Restrict periods to selected markers');
    [datatmp,headertmp,Idx2] = lab_include_markers(datatmp,headertmp,cfg.CONNECT.MARKER.markerinclude);
    if isempty(datatmp)
        disp('     Abort: calculation of connectivity not possible');
        return
    end
    Dindex = Dindex(1,Idx2);
    clearvars Idx2
end

% prepare skip first timeframes if selected
if isfield(cfg.CONNECT.PHASE,'skipstart') & cfg.CONNECT.PHASE.skipstart > 0 & cfg.CONNECT.PHASE.skipstart < (size(datatmp,2)-1)
    if size(datatmp,2) > cfg.CONNECT.PHASE.skipstart
        disp(['    Exclude first ' num2str(cfg.CONNECT.PHASE.skipstart) ' timeframes']);
        Dindex = Dindex((cfg.CONNECT.PHASE.skipstart + 1):end);
        datatmp = datatmp(:,(cfg.CONNECT.PHASE.skipstart + 1):end);
        if isfield(headertmp,'events') & isfield(headertmp.events,'POS') & ~isempty(headertmp.events.POS)
            headertmp.events.POS = headertmp.events.POS - cfg.CONNECT.PHASE.skipstart;
            tmp = find(headertmp.events.POS > 0);
            headertmp.events.POS = headertmp.events.POS(1,tmp);
            headertmp.events.DUR = headertmp.events.DUR(1,tmp);
            headertmp.events.OFF = headertmp.events.OFF(1,tmp);
            headertmp.events.TYP = headertmp.events.TYP(1,tmp);
            clearvars tmp
        end
        skipstart = cfg.CONNECT.PHASE.skipstart;
    else
        disp('     skip excluding first TFs, duration to short');
    end
else
    skipstart = 0;
end

% prepare skip last timeframes if selected
if isfield(cfg.CONNECT.PHASE,'skipend') & cfg.CONNECT.PHASE.skipend > 0 & cfg.CONNECT.PHASE.skipend < (size(datatmp,2)-1)
    if size(datatmp,2) > cfg.CONNECT.PHASE.skipend
        disp(['    Exclude last ' num2str(cfg.CONNECT.PHASE.skipend) ' timeframes']);
        Dindex = Dindex(1:(end-cfg.CONNECT.PHASE.skipend));
        datatmp = datatmp(1:(end-cfg.CONNECT.PHASE.skipend));
        skipend = cfg.CONNECT.PHASE.skipend;
        if isfield(headertmp,'events') & isfield(headertmp.events,'POS') & ~isempty(headertmp.events.POS)
            tmp = find(headertmp.events.POS <= size(datatmp,2));
            headertmp.events.POS = headertmp.events.POS(1,tmp);
            headertmp.events.DUR = headertmp.events.DUR(1,tmp);
            headertmp.events.OFF = headertmp.events.OFF(1,tmp);
            headertmp.events.TYP = headertmp.events.TYP(1,tmp);
            clearvars tmp
        end
    else
        disp('     skip excluding last TFs, duration to short');
    end
else
    skipend = 0;
end

% Select TFs by minimal correlation (only for data with microstates)
if isfield(headertmp,'CORR')
    headertmp.CORR = headertmp.CORR(:,Dindex);
    if isfield(cfg.CONNECT,'MARKER') & isfield(cfg.CONNECT.MARKER,'mincorr') & cfg.CONNECT.MARKER.mincorr > 0
        disp(['     restict to minimal cluster-correlation of ' num2str(cfg.CONNECT.MARKER.mincorr)]);
        tmp = abs(headertmp.CORR) >= cfg.CONNECT.MARKER.mincorr;
        if ~isempty(tmp)
            datatmp = datatmp(:,tmp);
            Dindex = Dindex(1,tmp);
            headertmp.events = lab_reduce_events(headertmp.events,tmp);
            headertmp.numtimeframes = size(datatmp,2);
            clearvars tmp
        else
            disp('     Abort: calculation of connectivity not possible');
            return
        end
    end
end

windowsize = cfg.CONNECT.PHASE.freqwindow;
if ~isfield(cfg.CONNECT.PHASE,'stepunit')
    cfg.CONNECT.PHASE.stepunit = 'TFs';
end
if cfg.CONNECT.PHASE.step > 0 & strcmp(cfg.CONNECT.PHASE.stepunit,'TFs') & cfg.CONNECT.PHASE.step <= size(datatmp,2);
    step = cfg.CONNECT.PHASE.step;
elseif cfg.CONNECT.PHASE.step > 0 & strcmp(cfg.CONNECT.PHASE.stepunit,'seconds') & ...
        isfield(headertmp,'samplingrate') & cfg.CONNECT.PHASE.step*headertmp.samplingrate <= size(datatmp,2);
    step = cfg.CONNECT.PHASE.step * headertmp.samplingrate;
elseif cfg.CONNECT.PHASE.step > 0 & strcmp(cfg.CONNECT.PHASE.stepunit,'phases of max frequency')
    step = floor(headertmp.samplingrate * (cfg.CONNECT.PHASE.step / max(Freqs(:))));
else
    step = size(datatmp,2);
end
if step < 60 % step should not drop below 60
    step = 60;
end
if step > size(datatmp,2);
    step = size(datatmp,2);
end
if windowsize < (step + skipstart + skipend)
    windowsize = step + skipstart + skipend;
end
if mod(windowsize,2) ~= 0
    windowsize = windowsize + 1;
end
stepmax = floor((size(datatmp,2) - (windowsize + step + skipstart)) / step);
if stepmax < 1
    stepmax = 1;
    step = size(datatmp,2);
end
skipstart = skipstart + floor((windowsize - (step+skipstart+skipend)) / 2);
skipend = skipend + ceil((windowsize - (step+skipstart+skipend)) / 2); %#ok<NASGU>

% find bad epochs
if isfield(cfg.CONNECT,'EPOCH') & ~isempty(cfg.CONNECT.EPOCH) & stepmax > 1
    disp('     Find epochs with automated routines');
    cfg.CONNECT.EPOCH.length = step;
    cfg.CONNECT.EPOCH = lab_markerselect(cfg.CONNECT.EPOCH,cfg,headertmp);
    [Epochs,cfg.CONNECT.EPOCH] = lab_select_epochs(datatmp,headertmp,cfg.CONNECT.EPOCH,windowsize);
    if ~isempty(Epochs)
        disp(['     ' num2str(size(Epochs.events.POS,2)) ' epochs found']);
        Dtmp = zeros(size(Epochs.events.POS,2),step);
        for Nepoch = 1:size(Epochs.events.POS,2)
            Dtmp(Nepoch,:) = Dindex(Epochs.events.POS(Nepoch):(Epochs.events.POS(Nepoch)+windowsize-1));
        end
        Dindex = Dtmp;
        clearvars Dtmp
        stepmax = size(Dindex,1);
    else
        disp('      find epochs failed, process all data');
        Dtmp = zeros(stepmax,windowsize);
        for nstep = 1:stepmax
            Dtmp(nstep,:) = Dindex(:,(nstep-1)*step + 1:(nstep-1)*step + windowsize);
        end
        Dindex = Dtmp;
        clearvars Dtmp
    end
else
    Dtmp = zeros(stepmax,windowsize);
    for nstep = 1:stepmax
        Dtmp(nstep,:) = Dindex(:,(nstep-1)*step + 1:(nstep-1)*step + windowsize);
    end
    Dindex = Dtmp;
    clearvars Dtmp
end

% Prepare Hanning-window
hanningwin = ones(1,windowsize);
if isfield(cfg.CONNECT.PHASE,'usehann') & cfg.CONNECT.PHASE.usehann == 1
    hanntmp = hann((ceil(windowsize/2)))';
    hannS = ceil(size(hanntmp,2) / 2);
    hanningwin(1,1:hannS) = hanntmp(1,1:hannS);
    hanningwin(1,end-hannS:end) = hanntmp(1,end-hannS:end);
    clearvars hannS hanntmp
end

% Define frequencies
F = -header.samplingrate/2:header.samplingrate/windowsize:header.samplingrate/2;
F = F(1:windowsize);
lowfreq = find(F >= min(min(Freqs)), 1);
highfreq = find(F < max(max(Freqs)), 1, 'last' );
F = F(1,lowfreq:highfreq);
for i = 1:size(Freqs,1)
    freqs(i,1) = find(F >= Freqs(i,1), 1); %#ok<AGROW>
    freqs(i,2) = find(F < Freqs(i,2), 1,'last'); %#ok<AGROW>
end
clearvars i

if isfield(cfg,'Output_filepath')
    Connect_filepath = cfg.Output_filepath;
end

% Calculate TFR
Result = [];
tfr = zeros(size(data,1),highfreq-lowfreq+1,step);
for nstep = 1:stepmax
    Dtmp = Dindex(nstep,:);
    disp (['   Step ' num2str(nstep) ' of ' num2str(stepmax)])
    
    if isfield(header,'timestamp') & isnumeric(header.timestamp)
        timestamp = {[cfg.Output_fileS '_' num2str(header.timestamp + Dtmp(1) - 1)]};
    elseif stepmax > 1 & isfield(cfg,'Output_fileS')
        timestamp = {[cfg.Output_fileS '_' num2str(Dtmp(1))]};
    elseif isfield(cfg,'Output_fileS')
        timestamp = {cfg.Output_fileS};
    else
        timestamp = {num2str(Dtmp(1))};
    end
    
    % Calculate TFR
    datastep = data(:,min(Dtmp):max(Dtmp));
    datastep(:,1:windowsize/2) = datastep(:,1:windowsize/2) .* repmat(hanningwin(1:windowsize/2),size(data,1),1);
    datastep(:,windowsize/2+1:end) = datastep(:,windowsize/2+1:end) .* repmat(hanningwin(windowsize/2+1:end),size(data,1),1);
    for j = 1:size(datastep,1)
        disp (['   Calculate rihaczek-distribution chan ' num2str(j)])
        tmp = lab_rihaczek(datastep(j,:)', header.samplingrate, windowsize);
        tmp = tmp(:,Dtmp - min(Dtmp) + 1);
        tfr(j,:,:) = tmp(lowfreq:highfreq,skipstart+1:skipstart+step);
        clear tmp
    end
    clear datastep
    
    % Calculate PLI/dPLI/wPLI/PLV/wPLV/PLT
    disp (['   Calculate matrix ' num2str(nstep)])
    for Nfreq = 1:size(freqs,1)
        Outputname = ['F' num2str(Freqs(Nfreq,3)) 'F' num2str(Freqs(Nfreq,4))];
        
        % set max phase error (= error of 1 timeframe shift)
        if ~isfield(cfg.CONNECT.PHASE,'maxerror') | isempty(cfg.CONNECT.PHASE.maxerror)
            cfg.CONNECT.PHASE.maxerror = 2;
        end
        maxphaseerror = Freqs(Nfreq,2) / header.samplingrate * pi * 2 * cfg.CONNECT.PHASE.maxerror;
        if maxphaseerror > pi/4;
            maxphaseerror = pi/4;
            disp('    value for phase error to large, set to pi/4');
        end
    
        tfrtmp = tfr(:,freqs(Nfreq,1):freqs(Nfreq,2),:);
        if Nfreq == size(freqs,1)
            clearvars tfr
        end
        nchans = size(tfrtmp,1);
        for i = 1:length(cfg.CONNECT.measure)
            measure = cfg.CONNECT.measure{i};
            switch measure
                case 'PLI'
                    pli = zeros(nchans,nchans);
                case 'dPLI'
                    dpli = zeros(nchans,nchans);
                case 'PLV'
                    plv = zeros(nchans,nchans);
                case 'PLT'
                    plt = zeros(nchans,nchans);
                case 'wPLI'
                    wpli = zeros(nchans,nchans);
                case 'wPLV'
                    wplv = zeros(nchans,nchans);
                    if isfield(cfg.CONNECT.PHASE,'storeangle') & cfg.CONNECT.PHASE.storeangle == true
                        wplv_angle = zeros(nchans,nchans);
                    end
                case 'PTE'
                    disp('    skip PTE, calculation for rihaczek-distribution not possible');
            end
        end
        for j = 1:nchans
            tmp = tfrtmp .* repmat(tfrtmp(j,:,:).^-1,[nchans 1 1]);
            freqdiff = sum(tmp,2);
            freqdiff = permute(freqdiff,[1 3 2]);
            tmp2 = ~isnan(freqdiff);
            tmp2 = reshape(tmp2,size(freqdiff));
            tmp2 = sum(tmp2,1);
            tmp2 = (tmp2(:)>0);
            freqdiff = freqdiff(:,tmp2);
            clearvars tmp2
            phasediff = angle(freqdiff);
            phasediff(abs(phasediff) < 10^-12) = 0; % not really needed
            for i = 1:length(cfg.CONNECT.measure)
                measure = cfg.CONNECT.measure{i};
                switch measure
                    case 'PLI'
                        pli(j,:) = abs(mean(sign(phasediff),2));
                        % pli(j,:) = abs(sum(sign(phasediff),2) ./ sum(abs(sign(phasediff)),2));
                        pli(j,j) = 0;
                    case 'dPLI'
                        dpli(j,:) = mean((sign(phasediff)+1),2) ./ 2;
                        % dpli(j,:) = sum((sign(phasediff)+1),2) ./ (2*sum(abs(sign(phasediff)),2));
                        dpli(j,j) = 0;
                    case 'wPLI'
                        if isfield(cfg.CONNECT.PHASE,'debiased') & cfg.CONNECT.PHASE.debiased == true
                            wpli(j,:) = (abs(sum(diffangle,2)).^2 - sum(diffangle.^2,2)) ./ (sum(abs(diffangle),2).^2 - sum(diffangle.^2,2));
                        else
                            wpli(j,:) = abs(sum(phasediff,2)) ./ sum(abs(phasediff),2);
                        end
                        wpli(j,j) = 0;
                        wpli(j,wpli(j,:)<0) = 0;
                    case 'PLV'
                        plv(j,:) = abs(sum((freqdiff ./ abs(freqdiff)),2)) ./ size(freqdiff,2);
                        plv(j,j) = 0;
                    case 'wPLV'
                        tmp = abs(sum((freqdiff ./ abs(freqdiff)),2)) ./ size(freqdiff,2);
                        tmp(abs(sin(angle(sum((freqdiff ./ abs(freqdiff)),2)))) < sin(maxphaseerror)) = 0;
                        wplv(j,:) = tmp;
                        wplv(j,j) = 0;
                        if isfield(cfg.CONNECT.PHASE,'storeangle') & cfg.CONNECT.PHASE.storeangle == true
                            wplv_angle(j,:) = abs(sin(angle(sum((freqdiff ./ abs(freqdiff)),2))));
                            wplv_angle(j,j) = 0;
                        end
                    case 'PLT'
                        phasediff = sign(sin(phasediff));
                        phasediff(phasediff == 0) = 1;
                        phasediff = sum(diff(phasediff,[],2)~=0,2);
                        plt(j,:) = 1 - exp(-(phasediff+1).^-1*(size(freqdiff,2)/header.samplingrate));
                        plt(j,j) = 0;
                end
            end
        end
        clearvars tfrtmp
        
        % Load existing data from file
        if isfield(cfg,'Output_file') & ~isempty(cfg.Output_file)
            if isfield(cfg,'patient') & ~isempty(cfg.patient)
                Connect_file = [cfg.patient '_Conn_' Outputname '.mat'];
            else
                Connect_file = ['Conn_' Outputname '.mat'];
                cfg.patient = [];
            end
            if isfield(cfg.CONNECT,'deleteold') & islogical(cfg.CONNECT.deleteold) & cfg.CONNECT.deleteold == true
                warning off %#ok<WNOFF>
                if ~isfield(cfg,'listold')
                    if exist(fullfile(Connect_filepath,Connect_file),'file')
                        disp('     delete Conn_F...mat from previous run')
                        delete(fullfile(Connect_filepath,Connect_file));
                    end
                    cfg.listold = cellstr(fullfile(Connect_filepath,Connect_file));
                elseif ~strcmp(cfg.listold,fullfile(Connect_filepath,Connect_file))
                    if exist(fullfile(Connect_filepath,Connect_file),'file')
                        disp('     delete Conn_F...mat from previous run')
                        delete(fullfile(Connect_filepath,Connect_file));
                    end
                    cfg.listold = [cfg.listold cellstr(fullfile(Connect_filepath,Connect_file))];
                end
                warning on %#ok<WNON>
            end
            if exist(fullfile(Connect_filepath,Connect_file),'file')
                load(fullfile(Connect_filepath,Connect_file));
            end
        end
        if ~exist('result','var') | ~isfield(result,'PHASE')
            result.pli = [];
            result.dpli = [];
            result.wpli = [];
            result.plv = [];
            result.wpli = [];
            result.plt = [];
            result.pli_timestamp = [];
            result.dpli_timestamp = [];
            result.wpli_timestamp = [];
            result.plv_timestamp = [];
            result.wplv_timestamp = [];
            result.wplv_angle = [];
            result.plt_timestamp = [];
            result.PHASE.step = [];
            result.PHASE.firststep = [];
            result.PHASE.laststep = [];
            result.PHASE.phaseintime = [];
            if isfield(header,'channels')
                result.channels = header.channels;
            end
            if isfield(header,'locs')
                result.locs = header.locs;
            end
        end
        
        % Collect results
        for i = 1:length(cfg.CONNECT.measure)
            measure = cfg.CONNECT.measure{i};
            switch measure
                case 'PLI'
                    result.pli(:,:,end+1) = pli;
                    result.pli_timestamp = [result.pli_timestamp timestamp];
                case 'dPLI'
                    result.dpli(:,:,end+1) = dpli;
                    result.dpli_timestamp = [result.dpli_timestamp timestamp];
                case 'wPLI'
                    result.wpli(:,:,end+1) = wpli;
                    result.wpli_timestamp = [result.wpli_timestamp timestamp];
                case 'PLV'
                    result.plv(:,:,end+1) = plv;
                    result.plv_timestamp = [result.plv_timestamp timestamp];
                case 'wPLV'
                    result.wplv(:,:,end+1) = wplv;
                    result.wplv_timestamp = [result.wplv_timestamp timestamp];
                    if isfield(cfg.CONNECT.PHASE,'storeangle') & cfg.CONNECT.PHASE.storeangle == true
                        result.wplv_angle(:,:,end+1) = wplv_angle;
                    end
                case 'PLT'
                    result.plt(:,:,end+1) = plt;
                    result.plt_timestamp = [result.plt_timestamp timestamp];
            end
        end
        result.freqband = Freqs(Nfreq,3:4);
        result.PHASE.freqs = Freqs(Nfreq,3:4);
        result.PHASE.step = step;
        result.PHASE.firststep = datstart;
        result.PHASE.laststep = (stepmax-1)*step + datstart;
        if isfield(header,'channels')
            result.channels = header.channels;
        end
        
        if exist('Connect_file','var')
            save(fullfile(Connect_filepath,Connect_file),'result');
        end
        
        if nstep == stemax
            Result.(Outputname) = result;
        end
        clearvars result
    end
end