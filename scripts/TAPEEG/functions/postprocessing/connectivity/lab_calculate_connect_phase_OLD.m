% Calculate PLI,dPLI,wPLI,num2str,wPLV,PLT and PTE based in specified frequency ranges
%
% [Result,cfg] = lab_calculate_connect_phase(data,header,cfg)
%
% data                = Matrix channels x timepoints
% header.samplingrate = Samplingrate
% cfg.CONNECT.step    = Size of steps in phases of max frequency
% cfg.CONNECT.freqs   = Matrix with frequencybands to analyze (low1 high1;low2 high2;...)
%
% written by F.Hatz 2012 Vumc Neurophysiology Amsterdam

function [Result,cfg,header] = lab_calculate_connect_phase(data,header,cfg)

if isfield(cfg,'Output_file') & ~isempty(cfg.Output_file)
    [~,~,~,cfg.Output_fileS] = lab_filename(cfg.Output_file);
end

if ~exist('cfg','var') || ~isfield(cfg,'CONNECT') || ~isfield(cfg.CONNECT,'freqs')
    cfg.CONNECT.freqs = [1 4;4 8;8 10;10 13;13 30];
end
if ~isfield(cfg.CONNECT,'PHASE') | ~isfield(cfg.CONNECT.PHASE,'step')
    cfg.CONNECT.PHASE.step = 8;  % phases of max freq
end
if ~isfield(cfg.CONNECT.PHASE,'window')
    cfg.CONNECT.PHASE.window = 4;  % min window for analysis in seconds
end
if ~isfield(cfg.CONNECT.PHASE,'freqwindow')
    if isfield(cfg.CONNECT.PHASE,'window') & cfg.CONNECT.PHASE.window == 0
        cfg.CONNECT.PHASE.freqwindow = size(data,2);
    else
        cfg.CONNECT.PHASE.freqwindow = floor(header.samplingrate*cfg.CONNECT.PHASE.window);
    end
end

if strcmp(cfg.CONNECT.filter,'input')
    if isempty(cfg.CONNECT.freqs) | size(cfg.CONNECT.freqs,2) ~= 2
        cfg.CONNECT.freqs = [0 0];
    else
        cfg.CONNECT.freqs = cfg.CONNECT.freqs(1,:);
    end
end

Result = [];
IdxAll = [];
Epochs = [];
for Nfreq = 1:size(cfg.CONNECT.freqs,1)
    headertmp = header;
    if max(cfg.CONNECT.freqs(Nfreq,:)) == 0
        disp ('   Calculate Connectivity unfiltered signal')
        Outputname = 'UNFILTERED';
    else
        disp (['   Calculate Connectivity ' num2str(cfg.CONNECT.freqs(Nfreq,1)) 'Hz to ' num2str(cfg.CONNECT.freqs(Nfreq,2)) 'Hz'])
        Outputname = ['F' num2str(cfg.CONNECT.freqs(Nfreq,1)) 'F' num2str(cfg.CONNECT.freqs(Nfreq,2))];
    end
    
    % filter data
    cfgfilt.lowpass = cfg.CONNECT.freqs(Nfreq,2);
    cfgfilt.highpass = cfg.CONNECT.freqs(Nfreq,1);
    if strcmp(cfg.CONNECT.filter,'freq-shift')
        disp('     set shifting for frequency filter')
        cfgfilt.filtermode = 'freq';
        cfgfilt.dofreqshift = true;
    else
        cfgfilt.filtermode = cfg.CONNECT.filter;
        cfgfilt.dofreqshift = false;
    end
    
    % set max phase error (= error of 1 timeframe shift)
    if ~isfield(cfg.CONNECT.PHASE,'maxerror') | isempty(cfg.CONNECT.PHASE.maxerror)
        cfg.CONNECT.PHASE.maxerror = 2;
    end
    maxphaseerror = cfg.CONNECT.freqs(Nfreq,2) / headertmp.samplingrate * pi * 2 * cfg.CONNECT.PHASE.maxerror;
    if maxphaseerror > pi/4;
        maxphaseerror = pi/4;
        disp('     value for phase error to large, set to pi/4');
    end
    
    % calculate phase
    if strcmp(cfg.CONNECT.PHASE.phaseestimate,'hilbert')
        phaseintime = lab_hilbert(data,headertmp,cfg.CONNECT.PHASE,cfgfilt);
        lag = 0;
    elseif strcmp(cfg.CONNECT.PHASE.phaseestimate,'zerocross')
        [phaseintime,lag] = lab_zerocross(data,headertmp,cfgfilt);
    elseif strcmp(cfg.CONNECT.PHASE.phaseestimate,'realphase')
        phaseintime = lab_calculate_realphase(data,headertmp,cfgfilt);
        lag = 0;
    end
    
    % correct header(events and CORR) for lag
    if lag > 0
        if isfield(headertmp,'events') & isfield(headertmp.events,'POS') & ~isempty(headertmp.events.POS)
            headertmp.events.POS = headertmp.events.POS - lag;
            tmp = intersect(find(headertmp.events.POS <= size(phaseintime,2)),find(headertmp.events.POS > 0));
            headertmp.events.POS = headertmp.events.POS(1,tmp);
            headertmp.events.DUR = headertmp.events.DUR(1,tmp);
            headertmp.events.OFF = headertmp.events.OFF(1,tmp);
            headertmp.events.TYP = headertmp.events.TYP(1,tmp);
            clearvars tmp
        end
        if isfield(headertmp,'CORR')
            headertmp.CORR = headertmp.CORR(:,lag+1:size(phaseintime,2));
        end
        headertmp.numtimeframes = size(phaseintime,2);
    end
    
    if isempty(IdxAll)
        % Create indexes
        Idx = 1:size(phaseintime,2);
        
        % exclude markers if necessary
        if isfield(cfg.CONNECT,'MARKER') & isfield(cfg.CONNECT.MARKER,'markerexclude') & ...
                iscell(cfg.CONNECT.MARKER.markerexclude) & ~isempty(cfg.CONNECT.MARKER.markerexclude) & ...
                isfield(headertmp,'events') & ~isempty(headertmp.events)
            disp('     exclude phase-periods with selected markers');
            [phaseintime,headertmp,Idx2] = lab_exclude_markers(phaseintime,headertmp,cfg.CONNECT.MARKER.markerexclude);
            if isempty(phaseintime)
                disp('     Abort: calculation of connectivity not possible');
                return
            end
            Idx = Idx(1,Idx2);
            clearvars Idx2
        end
        
        % include markers if necessary
        if isfield(cfg.CONNECT,'MARKER') & isfield(cfg.CONNECT.MARKER,'markerinclude') & ...
                iscell(cfg.CONNECT.MARKER.markerinclude) & ~isempty(cfg.CONNECT.MARKER.markerinclude) & ...
                isfield(headertmp,'events') & ~isempty(headertmp.events)
            disp('     restrict phase-periods to selected markers');
            [phaseintime,headertmp,Idx2] = lab_include_markers(phaseintime,headertmp,cfg.CONNECT.MARKER.markerinclude);
            if isempty(phaseintime)
                disp('     Abort: calculation of connectivity not possible');
                return
            end
            Idx = Idx(1,Idx2);
            clearvars Idx2
        end
        
        % skip first timeframes if selected
        if isfield(cfg.CONNECT.PHASE,'skipstart') & cfg.CONNECT.PHASE.skipstart > 0 & cfg.CONNECT.PHASE.skipstart < (size(phaseintime,2)-1)
            if size(phaseintime,2) > cfg.CONNECT.PHASE.skipstart
                disp(['     exclude first ' num2str(cfg.CONNECT.PHASE.skipstart) ' timeframes of phase information']);
                phaseintime = phaseintime(:,(cfg.CONNECT.PHASE.skipstart + 1):end);
                if isfield(headertmp,'events') & isfield(headertmp.events,'POS') & ~isempty(headertmp.events.POS)
                    headertmp.events.POS = headertmp.events.POS - cfg.CONNECT.PHASE.skipstart;
                    tmp = find(headertmp.events.POS > 0);
                    headertmp.events.POS = headertmp.events.POS(1,tmp);
                    headertmp.events.DUR = headertmp.events.DUR(1,tmp);
                    headertmp.events.OFF = headertmp.events.OFF(1,tmp);
                    headertmp.events.TYP = headertmp.events.TYP(1,tmp);
                    clearvars tmp
                end
                Idx = Idx(1,(cfg.CONNECT.PHASE.skipstart + 1):end);
            else
                disp('     skip excluding first TFs, duration to short');
            end
        end
        
        % skip last timeframes if selected
        if isfield(cfg.CONNECT.PHASE,'skipend') & cfg.CONNECT.PHASE.skipend > 0 & cfg.CONNECT.PHASE.skipend < (size(phaseintime,2)-1)
            if size(phaseintime,2) > cfg.CONNECT.PHASE.skipend
                disp(['     exclude last ' num2str(cfg.CONNECT.PHASE.skipend) ' timeframes of phase information']);
                phaseintime = phaseintime(:,1:(end-cfg.CONNECT.PHASE.skipend));
                if isfield(headertmp,'events') & isfield(headertmp.events,'POS') & ~isempty(headertmp.events.POS)
                    tmp = find(headertmp.events.POS <= size(phaseintime,2));
                    headertmp.events.POS = headertmp.events.POS(1,tmp);
                    headertmp.events.DUR = headertmp.events.DUR(1,tmp);
                    headertmp.events.OFF = headertmp.events.OFF(1,tmp);
                    headertmp.events.TYP = headertmp.events.TYP(1,tmp);
                    clearvars tmp
                end
                Idx = Idx(1,1:(end-cfg.CONNECT.PHASE.skipend));
            else
                disp('     skip excluding last TFs, duration to short');
            end
        end
        
        % Select TFs by minimal correlation (only for data with microstates)
        if isfield(headertmp,'CORR')
            headertmp.CORR = headertmp.CORR(:,Idx+lag);
            if isfield(cfg.CONNECT,'MARKER') & isfield(cfg.CONNECT.MARKER,'mincorr') & cfg.CONNECT.MARKER.mincorr > 0
                disp(['     restict to minimal cluster-correlation of ' num2str(cfg.CONNECT.MARKER.mincorr)]);
                Idx2 = abs(headertmp.CORR) >= cfg.CONNECT.MARKER.mincorr;
                if ~isempty(Idx2)
                    phaseintime = phaseintime(:,Idx2);
                    headertmp.CORR = headertmp.CORR(:,Idx2);
                    if isfield(headertmp,'events') & isfield(headertmp.events,'POS') & ~isempty(headertmp.events.POS)
                        headertmp.events = lab_reduce_events(headertmp.events,Idx2);
                    end
                    Idx = Idx(1,Idx2);
                    clearvars Idx2
                else
                    disp('     Abort: calculation of connectivity not possible');
                    return
                end
            end
        end
        IdxAll = Idx;
    else
        Idx = IdxAll;
        phaseintime = phaseintime(:,Idx);
        if isfield(headertmp,'CORR')
            headertmp.CORR = headertmp.CORR(:,Idx);
        end
        if isfield(headertmp,'events') & isfield(headertmp.events,'POS') & ~isempty(headertmp.events.POS)
            headertmp.events = lab_reduce_events(headertmp.events,Idx);
        end
    end
    
    % set step
    if ~isfield(cfg.CONNECT.PHASE,'stepunit')
        cfg.CONNECT.PHASE.stepunit = 'TFs';
    end
    if cfg.CONNECT.PHASE.step > 0 & strcmp(cfg.CONNECT.PHASE.stepunit,'TFs') & cfg.CONNECT.PHASE.step <= size(phaseintime,2);
        step = cfg.CONNECT.PHASE.step;
    elseif cfg.CONNECT.PHASE.step > 0 & strcmp(cfg.CONNECT.PHASE.stepunit,'seconds') & ...
            isfield(headertmp,'samplingrate') & cfg.CONNECT.PHASE.step*headertmp.samplingrate <= size(phaseintime,2);
        step = cfg.CONNECT.PHASE.step * headertmp.samplingrate;
    elseif cfg.CONNECT.PHASE.step > 0 & strcmp(cfg.CONNECT.PHASE.stepunit,'phases of max frequency')
        step = floor(headertmp.samplingrate * (cfg.CONNECT.PHASE.step / max(cfg.CONNECT.freqs(Nfreq,:))));
    else
        step = size(phaseintime,2);
    end
    if step < 60 % step should not drop below 60
        step = 60;
    end
    if step > size(phaseintime,2);
        step = size(phaseintime,2);
    end
    stepmax = floor(size(phaseintime,2) / step);
    if stepmax < 1
        stepmax = 1;
        step = size(phaseintime,2);
    end
    disp(['     step is set to ' num2str(step) ', length is ' num2str(size(phaseintime,2))]);
    
    % find bad epochs
    if isfield(cfg.CONNECT,'EPOCH') & ~isempty(cfg.CONNECT.EPOCH) & stepmax > 1
        if isempty(Epochs) | size(Epochs.data,2) ~= step
            disp('     Find epochs with automated routines');
            headertmp.numtimeframes = length(Idx);
            cfg.CONNECT.EPOCH.length = step;
            cfg.CONNECT.EPOCH = lab_markerselect(cfg.CONNECT.EPOCH,cfg,headertmp);
            [Epochs,cfg.CONNECT.EPOCH] = lab_select_epochs(data(:,Idx+lag),headertmp,cfg.CONNECT.EPOCH);
        end
        if ~isempty(Epochs)
            disp(['     ' num2str(size(Epochs.events.POS,2)) ' epochs found']);
            phasetmp = zeros(size(phaseintime,1),step,size(Epochs.events.POS,2));
            for Nepoch = 1:size(Epochs.events.POS,2)
                phasetmp(:,:,Nepoch) = phaseintime(:,Epochs.events.POS(Nepoch):(Epochs.events.POS(Nepoch)+step-1));
            end
            phaseintime = phasetmp;
            clearvars phasetmp
            stepmax = size(phaseintime,3);
        else
            disp('      find epochs failed, process all data');
            phaseintime = phaseintime(:,1:stepmax*step);
            phaseintime = reshape(phaseintime,[size(phaseintime,1) step stepmax]);
        end
    else
        phaseintime = phaseintime(:,1:stepmax*step);
        phaseintime = reshape(phaseintime,[size(phaseintime,1) step stepmax]);
        if isfield(headertmp,'CORR')
            CORR = headertmp.CORR(1,1:stepmax*step);
            CORR = reshape(CORR,[1 step stepmax]);
            CORR = mean(permute(abs(CORR),[3 2 1]),2);
            [~,Cidx] = sort(CORR,'descend');
            phaseintime = phaseintime(:,:,Cidx);
            clearvars CORR Cidx
        end
    end
    
    % prepare calculation
    result = [];
    DispMeasure = '';
    Nchans = size(phaseintime,1);
    for i = 1:length(cfg.CONNECT.measure)
        measure = cfg.CONNECT.measure{i};
        switch measure
            case 'PLI'
                result.pli = zeros(Nchans,Nchans,stepmax);
                DispMeasure = [DispMeasure 'PLI ']; %#ok<AGROW>
            case 'dPLI'
                result.dpli = zeros(Nchans,Nchans,stepmax);
                DispMeasure = [DispMeasure 'dPLI ']; %#ok<AGROW>
            case 'wPLI'
                result.wpli = zeros(Nchans,Nchans,stepmax);
                DispMeasure = [DispMeasure 'wPLI ']; %#ok<AGROW>
            case 'PLV'
                result.plv = zeros(Nchans,Nchans,stepmax);
                DispMeasure = [DispMeasure 'PLV ']; %#ok<AGROW>
            case 'wPLV'
                result.wplv = zeros(Nchans,Nchans,stepmax);
                if isfield(cfg.CONNECT.PHASE,'storeangle') & cfg.CONNECT.PHASE.storeangle == true
                    result.wplv_angle = zeros(Nchans,Nchans,stepmax);
                end
                DispMeasure = [DispMeasure 'wPLV ']; %#ok<AGROW>
            case 'PLT'
                result.plt = zeros(Nchans,Nchans,stepmax);
                DispMeasure = [DispMeasure 'PLT ']; %#ok<AGROW>
            case 'PTE'
                result.pte = zeros(Nchans,Nchans,stepmax);
                pte_ci = zeros(Nchans,Nchans,stepmax);
                DispMeasure = [DispMeasure 'PTE ']; %#ok<AGROW>
                if ~isfield(cfg.CONNECT.PHASE,'delta') | isempty(cfg.CONNECT.PHASE.delta)
                    disp('     PTE: factor for delta not set, set to default = 2')
                    cfg.CONNECT.PHASE.delta = 2;
                end
                H = zeros(Nchans,1);
                % D = [];
                for j = 1:Nchans
                    tmp = phaseintime(j,:,:);
                    H(j,1) = (3.5 * circ_std(tmp(:))) / (length(tmp(:))^(1/3));
                    % D = cat(2,D,diff(find(abs(diff(sign(tmp)))==2)));
                end
                clearvars tmp
                H = median(H);
                H = round((2*pi/H)/2)*2;
                phaseintimeH = round((phaseintime + pi)/(2*pi) * (H - 1)) + 1;
                
                if isfield(cfg.CONNECT.PHASE,'deltafreq') & cfg.CONNECT.PHASE.deltafreq == true
                    delta = round((header.samplingrate / ((cfgfilt.lowpass + cfgfilt.highpass)/2)) / cfg.CONNECT.PHASE.delta);
                else
                    delta = cfg.CONNECT.PHASE.delta;
                end
                % delta = round(mean(D) / cfg.CONNECT.PHASE.delta);
                if delta < 1
                    delta = 1;
                end
                disp(['     PTE delta = ', num2str(delta)])
        end
    end
    
    % calculate connectivity
    
    for nstep = 1:stepmax
        disp(['     calculate ' DispMeasure '(matrix ' num2str(nstep) ' of ' num2str(stepmax) ')'])
        phasetmp = phaseintime(:,:,nstep);
        phasetmp = ones(size(phasetmp,1),size(phasetmp,2)).*exp(1i*phasetmp);
        for nchan = 1:Nchans
            phasediff = phasetmp .* repmat(phasetmp(nchan,:).^-1,Nchans,1);
            diffangle = angle(phasediff);
            % diffangle(abs(diffangle) < 10^-12) = 0; % not really needed, but doesn't alter result
            for i = 1:length(cfg.CONNECT.measure)
                measure = cfg.CONNECT.measure{i};
                switch measure
                    case 'PLI'
                        result.pli(:,nchan,nstep) = abs(mean(sign(sin(diffangle)),2));
                        result.pli(nchan,nchan,nstep) = 0;
                    case 'dPLI'
                        result.dpli(:,nchan,nstep) = mean((sign(sin(diffangle))+1),2) ./ 2;
                        result.dpli(nchan,nchan,nstep) = 0;
                    case 'wPLI'
                        if isfield(cfg.CONNECT.PHASE,'debiased') & cfg.CONNECT.PHASE.debiased == true
                            result.wpli(:,nchan,nstep) = (abs(sum(diffangle,2)).^2 - sum(diffangle.^2,2)) ./ (sum(abs(diffangle),2).^2 - sum(diffangle.^2,2));
                        else
                            result.wpli(:,nchan,nstep) = abs(sum(diffangle,2)) ./ sum(abs(diffangle),2);
                        end
                        result.wpli(nchan,nchan,nstep) = 0;
                        result.wpli(result.wpli(:,nchan,nstep)<0,nchan,nstep) = 0;
                    case 'PLV'
                        result.plv(:,nchan,nstep) =  abs(sum((phasediff ./ abs(phasediff)),2)) ./ size(phasediff,2);
                        result.plv(nchan,nchan,nstep) = 0;
                    case 'wPLV'
                        tmp = (phasediff ./ abs(phasediff));
                        tmp(abs(sin(angle(tmp))) < sin(maxphaseerror)) = 0;
                        tmp = abs(sum(tmp,2)) ./ size(phasediff,2);
                        result.wplv(:,nchan,nstep) = tmp;
                        result.wplv(nchan,nchan,nstep) = 0;
                        clearvars tmp
                        if isfield(cfg.CONNECT.PHASE,'storeangle') & cfg.CONNECT.PHASE.storeangle == true
                            result.wplv_angle(:,nchan,nstep) = abs(sin(angle(sum((phasediff ./ abs(phasediff)),2))));
                            result.wplv_angle(nchan,nchan,nstep) = 0;
                        end
                    case 'PLT'
                        diffangle = sign(sin(diffangle));
                        diffangle(diffangle == 0) = 1;
                        for nchan2 = 1:Nchans
                            result.plt_shifts{nchan,nchan2,nstep} = find(diff(diffangle(nchan2,:),[],2)~=0);
                        end
                        result.plt_shifts{nchan,nchan,nstep} = find(diff(diffangle,[],2)~=0);
                        diffangle = sum(diff(diffangle,[],2)~=0,2);
                        result.plt(:,nchan,nstep) = 1 - exp(-(diffangle).^-1*(size(phasetmp,2)/headertmp.samplingrate));
                        result.plt(nchan,nchan,nstep) = 0;
                    case 'PTE'
                        phase1 = phaseintimeH(:,1+delta:end,nstep);
                        phase2 = phaseintimeH(:,1:end-delta,nstep);
                        L = size(phase2,2);
                        phaseT = (phase1-1)*H + phase2;
                        P1 = histc(phaseT(nchan,:),1:H^2,2) / L;
                        R1 = - nansum(P1 .* log(P1),2);
                        P3 = histc(phase2(nchan,:),1:H,2) / L;
                        R3 = - nansum(P3 .* log(P3),2);
                        P2 = histc((phase2-1)*H + repmat(phase2(nchan,:),Nchans,1),1:H^2,2) / L;
                        R2 = - nansum(P2 .* log(P2),2);
                        
                        P4 = histc((phaseT-1)*H + repmat(phase2(nchan,:),Nchans,1),1:H^3,2) / L;
                        R4 = - nansum(P4 .* log(P4),2);
                        result.pte(:,nchan,nstep) = R1 + R2 - R3 - R4;
                        result.pte(nchan,nchan,nstep) = 0;
                    case 'nPTE'
                        phase1 = phaseintimeH(:,1+delta:end,nstep);
                        phase2 = phaseintimeH(:,1:end-delta,nstep);
                        L = size(phase2,2);
                        phaseT = (phase1-1)*H + phase2;
                        P1 = histc(phaseT(nchan,:),1:H^2,2) / L;
                        R1 = - nansum(P1 .* log(P1),2);
                        P2 = histc((phase2-1)*H + repmat(phase2(nchan,:),Nchans,1),1:H^2,2) / L;
                        R2 = - nansum(P2 .* log(P2),2);
                        P3 = histc(phase2(nchan,:),1:H,2) / L;
                        R3 = - nansum(P3 .* log(P3),2); 
                        P4 = histc((phaseT-1)*H + repmat(phase2(nchan,:),Nchans,1),1:H^3,2) / L;
                        R4 = - nansum(P4 .* log(P4),2);
                        P5 = histc((phase1-1)*H + repmat(phase1(nchan,:),Nchans,1),1:H^2,2) / L;
                        R5 = - nansum(P5 .* log(P5),2);
                        result.pte(:,nchan,nstep) = R1 + R2 - R3 - R4;
                        result.pte(nchan,nchan,nstep) = 0;
                        pte_ic(:,nchan,nstep) = R5;
                end
            end
        end
        clearvars phasediff diffangle phasetmp nchan P1 P2 P3 P4 R1 R2 R3 R4 P5 R5 P6 R6 M phase1 phase2 L phaseT
        
        if isfield(headertmp,'timestamp') & isfield(cfg,'Output_fileS')
            timestamp{nstep} = [cfg.Output_fileS '_' num2str(headertmp.timestamp + (nstep-1)*step)]; %#ok<AGROW>
        elseif stepmax > 1 & isfield(cfg,'Output_fileS')
            timestamp{nstep} = [cfg.Output_fileS '_' num2str((nstep-1)*step + 1)]; %#ok<AGROW>
        elseif isfield(cfg,'Output_fileS')
            timestamp{nstep} = cfg.Output_fileS; %#ok<AGROW>
        else
            timestamp{nstep} = num2str((nstep-1)*step + 1); %#ok<AGROW>
        end
    end
    result.freqband = cfg.CONNECT.freqs(Nfreq,:);
    result.PHASE.freqs = cfg.CONNECT.freqs(Nfreq,:);
    result.PHASE.step = step;
    result.PHASE.firststep = lag + 1;
    result.PHASE.laststep = (stepmax-1)*step + lag + 1;
    if isfield(cfg.CONNECT.PHASE,'storephase') & cfg.CONNECT.PHASE.storephase == true
        result.PHASE.phaseintime = phaseintime;
    else
        result.PHASE.phaseintime = [];
    end
    if isfield(headertmp,'channels')
        result.PHASE.channels = headertmp.channels;
    end
    if isfield(result,'pli')
        result.pli_timestamp = timestamp;
    end
    if isfield(result,'dpli')
        result.dpli_timestamp = timestamp;
    end
    if isfield(result,'wpli')
        result.wpli_timestamp = timestamp;
    end
    if isfield(result,'plv')
        result.plv_timestamp = timestamp;
    end
    if isfield(result,'wplv')
        result.wplv_timestamp = timestamp;
    end
    if isfield(result,'plt')
        result.plt_timestamp = timestamp;
    end
    if isfield(result,'pte')
        result.pte_timestamp = timestamp;
        if isfield(cfg.CONNECT.PHASE,'donorm') & cfg.CONNECT.PHASE.donorm == true
            for i = 1:size(result.pte,3)
                tmp = result.pte(:,:,i);
                tmp = (tmp-tmp') ./ (tmp + tmp' + pte_ci(:,:,i));
                result.pte(:,:,i) = tmp;
            end
            clearvars tmp i;
        elseif isfield(cfg.CONNECT.PHASE,'dodirect') & cfg.CONNECT.PHASE.dodirect == true
            for i = 1:size(result.pte,3)
                tmp = result.pte(:,:,i);
                tmp = tmp ./ (tmp + tmp');
                result.pte(:,:,i) = tmp;
            end
            clearvars tmp i;
        end
    end
    if isfield(cfg,'Output_filepath') & ~isempty(cfg.Output_filepath)
        % save Result
        Connect_filepath = cfg.Output_filepath;
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
            elseif ~max(strcmp(cfg.listold,fullfile(Connect_filepath,Connect_file)))
                if exist(fullfile(Connect_filepath,Connect_file),'file')
                    disp('     delete Conn_F...mat from previous run')
                    delete(fullfile(Connect_filepath,Connect_file));
                end
                cfg.listold = [cfg.listold cellstr(fullfile(Connect_filepath,Connect_file))];
            end
            warning on %#ok<WNON>
        end
        if exist(fullfile(Connect_filepath,Connect_file),'file')
            result2 = result;
            clearvars result
            try %#ok<TRYNC>
                load(fullfile(Connect_filepath,Connect_file));
            end
            if ~exist('result','var') | ~isfield(result,'PHASE')
                result.pli = [];
                result.dpli = [];
                result.wpli = [];
                result.plv = [];
                result.wplv = [];
                result.plt = [];
                result.pte = [];
                result.pli_timestamp = [];
                result.dpli_timestamp = [];
                result.wpli_timestamp = [];
                result.plv_timestamp = [];
                result.wplv_timestamp = [];
                result.wplv_angle = [];
                result.plt_timestamp = [];
                result.pte_timestamp = [];
                result.plt_shift = [];
                result.wplv_angle = [];
                result.PHASE.step = [];
                result.PHASE.firststep = [];
                result.PHASE.laststep = [];
                result.PHASE.phaseintime = [];
            end
            if isfield(result2,'pli')
                result.pli = cat(3,result.pli,result2.pli);
                result.pli_timestamp = [result.pli_timestamp timestamp];
            end
            if isfield(result2,'dpli')
                result.dpli = cat(3,result.dpli,result2.dpli);
                result.dpli_timestamp = [result.dpli_timestamp timestamp];
            end
            if isfield(result2,'wpli')
                result.wpli = cat(3,result.wpli,result2.wpli);
                result.wpli_timestamp = [result.wpli_timestamp timestamp];
            end
            if isfield(result2,'plv')
                result.plv = cat(3,result.plv,result2.plv);
                result.plv_timestamp = [result.plv_timestamp timestamp];
            end
            if isfield(result2,'wplv')
                result.wplv = cat(3,result.wplv,result2.wplv);
                result.wplv_timestamp = [result.wplv_timestamp timestamp];
                if isfield(cfg.CONNECT.PHASE,'storeangle') & cfg.CONNECT.PHASE.storeangle == true
                    result.wplv_angle = cat(3,result.wplv_angle,result2.wplv_angle);
                end
            end
            if isfield(result2,'plt')
                result.plt = cat(3,result.plt,result2.plt);
                result.plt_shift = cat(2,result.plt_shift,result2.plt_shift);
                result.plt_timestamp = [result.plt_timestamp timestamp];
            end
            if isfield(result2,'pte')
                result.pte = cat(3,result.pte,result2.pte);
                result.pte_timestamp = [result.pte_timestamp timestamp];
            end
            result.PHASE.firststep = cat(1,result.PHASE.firststep,result2.PHASE.firststep);
            result.PHASE.laststep = cat(1,result.PHASE.laststep,result2.PHASE.laststep);
            result.PHASE.step = cat(1,result.PHASE.step,result2.PHASE.step);
            result.PHASE.phaseintime = [result.PHASE.phaseintime result2.PHASE.phaseintime];
        end
        if isfield(cfg,'patient')
            patient = cfg.patient; %#ok<NASGU>
        else
            patient = []; %#ok<NASGU>
        end
        if isfield(headertmp,'channels')
            result.channels = headertmp.channels;
        end
        if isfield(headertmp,'locs')
            result.locs = headertmp.locs;
        end
        result.freqband = cfg.CONNECT.freqs(Nfreq,:);
        save(fullfile(Connect_filepath,Connect_file),'result','patient');
    end
    Result.(Outputname) = result;
end