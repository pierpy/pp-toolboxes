% Calculate connectivity based on Morlet-wavelet analysis in
% specified frequency ranges
%
% [result,cfg] = lab_calculate_connect_wavelet(data,header,cfg)
%
% data                = Matrix channels x timepoints
% header.samplingrate = Samplingrate
% cfg.CONNECT.step    = Size of steps in phases of max frequency
% cfg.CONNECT.freqs   = Matrix with frequencybands to analyze (low1 high1;low2 high2;...)
% cfg.CONNECT.window  = Min window for wavelet analysis in seconds
%
% written by F.Hatz 2012 Vumc Neurophysiology Amsterdam
% (lab_specest_wavelet(..) from fieldtrip-toolbox)

function [Result,cfg,header] = lab_calculate_connect_wavelet(data,header,cfg)

if isfield(cfg,'Output_file')
    [~,~,~,cfg.Output_fileS] = lab_filename(cfg.Output_file);
end

if ~exist('cfg','var') || ~isfield(cfg,'CONNECT') || ~isfield(cfg.CONNECT,'freqs')
    cfg.CONNECT.freqs = [1 4;4 8;8 10;10 13;13 30];
end
if size(cfg.CONNECT.freqs,2) == 4
    Freqs = cfg.CONNECT.freqs;
elseif size(cfg.CONNECT.freqs,2) == 2
    Freqs = [cfg.CONNECT.freqs cfg.CONNECT.freqs];
else
    Freqs = [0 0 0 0];
end

if ~isfield(cfg.CONNECT,'PHASE') | ~isfield(cfg.CONNECT.PHASE,'step')
    cfg.CONNECT.PHASE.step = 8;  % phases of max freq
end
if ~isfield(cfg.CONNECT.PHASE,'window')
    cfg.CONNECT.PHASE.window = 4;  % min window for wavelet analysis in seconds
end
if ~isfield(cfg.CONNECT.PHASE,'freqwindow')
    cfg.CONNECT.PHASE.freqwindow = cfg.CONNECT.PHASE.window * header.samplingrate;  % window for wavelet analysis
end

Result = [];
IdxAll = [];
Epochs = [];
for Nfreq = 1:size(Freqs,1)
    disp (['   Calculate wavelets ' num2str(Freqs(Nfreq,3)) 'Hz to ' num2str(Freqs(Nfreq,4)) 'Hz'])
    
    
    % set temporary variables
    datatmp = data;
    headertmp = header;
    lag = ceil(header.samplingrate * 4 / min(Freqs(Nfreq,1:2)));
    
    if isempty(IdxAll)
        Dindex = 1:size(data,2);
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
            if cfg.CONNECT.PHASE.skipstart > lag
                skipstart = cfg.CONNECT.PHASE.skipstart;
            else
                skipstart = lag;
            end
        else
            skipstart = lag;
        end
        if skipstart < size(datatmp,2)
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
        else
            disp('     Abort: calculation of connectivity not possible');
            break
        end
        
        % prepare skip last timeframes
        if isfield(cfg.CONNECT.PHASE,'skipend') & cfg.CONNECT.PHASE.skipend > 0 & cfg.CONNECT.PHASE.skipend < (size(datatmp,2)-1)
            if cfg.CONNECT.PHASE.skipend > lag
                skipend = cfg.CONNECT.PHASE.skipend;
            else
                skipend = lag;
            end
        else
            skipend = lag;
        end
        if skipend < size(datatmp,2)
            disp(['    Exclude last ' num2str(skipend) ' timeframes']);
            Dindex = Dindex(1:(end-skipend));
            datatmp = datatmp(1:(end-skipend));
            if isfield(headertmp,'events') & isfield(headertmp.events,'POS') & ~isempty(headertmp.events.POS)
                tmp = find(headertmp.events.POS <= size(datatmp,2));
                headertmp.events.POS = headertmp.events.POS(1,tmp);
                headertmp.events.DUR = headertmp.events.DUR(1,tmp);
                headertmp.events.OFF = headertmp.events.OFF(1,tmp);
                headertmp.events.TYP = headertmp.events.TYP(1,tmp);
                clearvars tmp
            end
        else
            disp('     Abort: calculation of connectivity not possible');
            break
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
        IdxAll = Dindex;
    else
        Dindex = IdxAll;
        datatmp =datatmp(:,Dindex);
        if isfield(headertmp,'CORR')
            headertmp.CORR = headertmp.CORR(:,Dindex);
        end
        if isfield(headertmp,'events') & isfield(headertmp.events,'POS') & ~isempty(headertmp.events.POS)
            headertmp.events = lab_reduce_events(headertmp.events,Dindex);
        end
    end
    
    if cfg.CONNECT.PHASE.step > 0
        step = floor(header.samplingrate * (cfg.CONNECT.PHASE.step / max(Freqs(Nfreq,1:2))));
    else
        step = cfg.CONNECT.PHASE.freqwindow;
    end
    if step < 60 % step should not drop below 60
        step = 60;
    end
    if isempty(step)
        step = 100;
    end
    wavwindow = cfg.CONNECT.PHASE.window * header.samplingrate;
    if (2 * lag + step + skipstart + skipend) > wavwindow
        wavwindow = step + 2*lag + skipstart + skipend;
    end
    if mod(wavwindow,2) ~= 0
        wavwindow = wavwindow + 1;
    end
    stepmax = floor((size(datatmp,2) - (wavwindow + step + lag + skipstart)) / step);
    if stepmax < 1
        stepmax = 1;
        step = size(datatmp,2);
    end
    
    % find bad epochs
    if isfield(cfg.CONNECT,'EPOCH') & ~isempty(cfg.CONNECT.EPOCH) & stepmax > 1
        if isempty(Epochs) | size(Epochs.data,2) ~= step
            disp('     Find epochs with automated routines');
            cfg.CONNECT.EPOCH.length = step;
            cfg.CONNECT.EPOCH = lab_markerselect(cfg.CONNECT.EPOCH,cfg,headertmp);
            [Epochs,cfg.CONNECT.EPOCH] = lab_select_epochs(datatmp,headertmp,cfg.CONNECT.EPOCH,wavwindow);
        end
        if ~isempty(Epochs)
            disp(['     ' num2str(size(Epochs.events.POS,2)) ' epochs found']);
            Dtmp = zeros(size(Epochs.events.POS,2),step);
            for Nepoch = 1:size(Epochs.events.POS,2)
                Dtmp(Nepoch,:) = Dindex(Epochs.events.POS(Nepoch):(Epochs.events.POS(Nepoch)+wavwindow-1));
            end
            Dindex = Dtmp;
            clearvars Dtmp
            stepmax = size(Dindex,1);
        else
            disp('      find epochs failed, process all data');
            Dtmp = zeros(stepmax,wavwindow);
            for nstep = 1:stepmax
                Dtmp(nstep,:) = Dindex(:,(nstep-1)*step + 1:(nstep-1)*step + wavwindow);
            end
            Dindex = Dtmp;
            clearvars Dtmp
        end
    else
        Dtmp = zeros(stepmax,wavwindow);
        for nstep = 1:stepmax
            Dtmp(nstep,:) = Dindex(:,(nstep-1)*step + 1:(nstep-1)*step + wavwindow);
        end
        Dindex = Dtmp;
        clearvars Dtmp
    end
    
    % Define frequencies
    F = Freqs(Nfreq,1):header.samplingrate/wavwindow:Freqs(Nfreq,2);
    
    % set max phase error (= error of 1 timeframe shift)
    if ~isfield(cfg.CONNECT.PHASE,'maxerror') | isempty(cfg.CONNECT.PHASE.maxerror)
        cfg.CONNECT.PHASE.maxerror = 2;
    end
    maxphaseerror = Freqs(Nfreq,2) / header.samplingrate * pi * 2 * cfg.CONNECT.PHASE.maxerror;
    if maxphaseerror > pi/4;
        maxphaseerror = pi/4;
        disp('    value for phase error to large, set to pi/4');
    end
    
    % Calculate Wavelet and (d)PLI
    Outputname = ['F' num2str(Freqs(Nfreq,3)) 'F' num2str(Freqs(Nfreq,4))];
    time = 0:1/header.samplingrate:wavwindow/header.samplingrate;
    time = time(1,1:wavwindow);
    result.phaseintime = [];
    
    if isfield(cfg,'Output_filepath')
        % try load results
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
            elseif ~strcmp(cfg.listold,fullfile(Connect_filepath,Connect_file))
                if exist(fullfile(Connect_filepath,Connect_file),'file')
                    disp('     delete Conn_F...mat from previous run')
                    delete(fullfile(Connect_filepath,Connect_file));
                end
                cfg.listold = [cfg.listold cellstr(fullfile(Connect_filepath,Connect_file))];
            end
            warning on %#ok<WNON>
        end
        try %#ok<TRYNC>
            load(fullfile(Connect_filepath,Connect_file));
        end
    end
    if ~exist('result','var') | ~isfield(result,'PHASE')
        result.pli = [];
        result.dpli = [];
        result.wpli = [];
        result.plv = [];
        result.wplv = [];
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
    
    for nstep = 1:stepmax
        Dtmp = Dindex(nstep,:);
        disp (['   Step ' num2str(nstep) ' of ' num2str(stepmax)])
        
        if isfield(header,'timestamp')
            timestamp{nstep} = [cfg.Output_fileS '_' num2str(header.timestamp + Dtmp(1) - 1)]; %#ok<AGROW>
        elseif stepmax > 1 & isfield(cfg,'Output_fileS')
            timestamp{nstep} = [cfg.Output_fileS '_' num2str(Dtmp(1))]; %#ok<AGROW>
        elseif isfield(cfg,'Output_fileS')
            timestamp{nstep} = cfg.Output_fileS; %#ok<AGROW>
        else
            timestamp{nstep} = num2str(Dtmp(1)); %#ok<AGROW>
        end
        
        % Calculate Wavelet
        datastep = data(:,min(Dtmp):max(Dtmp));
        wavelettmp = lab_specest_wavelet(datastep,time,'freqoi',F);
        wavelettmp = wavelettmp(:,:,Dtmp - min(Dtmp) + 1);
        wavelettmp = wavelettmp(:,:,skipstart+1:skipstart+step);
        freqdiff = zeros(size(wavelettmp,1),size(wavelettmp,1),size(wavelettmp,3));
        for i = 1:size(wavelettmp,3)
            tmp = zeros(size(wavelettmp,1),size(wavelettmp,1),size(wavelettmp,2));
            for m = 1:size(wavelettmp,2)
                tmp(:,:,m) = wavelettmp(:,m,i) * ((wavelettmp(:,m,i).').^-1);
            end
            clearvars m
            freqdiff(:,:,i) = sum(tmp,3);
        end
        clearvars i wavelettmp
        tmp = ~isnan(freqdiff);
        tmp = reshape(tmp,size(freqdiff));
        tmp = sum(sum(tmp,1),2);
        tmp = (tmp(:)>0);
        freqdiff = freqdiff(:,:,tmp);
        clearvars tmp
        phasediff = angle(freqdiff);
        % phasediff(abs(phasediff) < 10^-12) = 0; % not really needed
        
        % Calculate PLI/dPLI/wPLI/PLV/wPLV/PLT
        for i = 1:length(cfg.CONNECT.measure)
            measure = cfg.CONNECT.measure{i};
            switch measure
                case 'PLI'
                    tmp = abs(mean(sign(phasediff),3));
                    % tmp = abs(sum(sign(phasediff),3) ./ sum(abs(sign(phasediff)),3));
                    tmp(1:size(tmp,1)+1:end) = 0;
                        result.pli(:,:,end+1) = tmp;
                        result.pli_timestamp = [result.pli_timestamp timestamp{nstep}];
                case 'dPLI'
                    tmp = mean((sign(phasediff)+1),3) ./ 2;
                    % tmp = sum((sign(phasediff)+1),3) ./ (2*sum(abs(sign(phasediff)),3));
                    tmp(1:size(tmp,1)+1:end) = 0;
                    result.dpli(:,:,end+1) = tmp;
                    result.dpli_timestamp = [result.dpli_timestamp timestamp{nstep}];
                case 'wPLI'
                    if isfield(cfg.CONNECT.PHASE,'debiased') & cfg.CONNECT.PHASE.debiased == true
                        tmp = (abs(sum(diffangle,3)).^2 - sum(diffangle.^2,3)) ./ (sum(abs(diffangle),3).^2 - sum(diffangle.^2,3));
                    else
                        tmp = abs(sum(phasediff,3)) ./ sum(abs(phasediff),3);
                    end
                    tmp(1:size(tmp,1)+1:end) = 0;
                    tmp(tmp < 0) = 0;
                    result.wpli(:,:,end+1) = tmp;
                    result.wpli_timestamp = [result.wpli_timestamp timestamp{nstep}];
                    clearvars tmp
                case 'PLV'
                    tmp = abs(sum((freqdiff ./ abs(freqdiff)),3)) ./ size(freqdiff,3);
                    tmp(1:size(tmp,1)+1:end) = 0;
                    result.plv(:,:,end+1) = tmp;
                    result.plv_timestamp = [result.plv_timestamp timestamp{nstep}];
                case 'wPLV'
                    tmp = abs(sum((freqdiff ./ abs(freqdiff)),3)) ./ size(freqdiff,3);
                    tmp(abs(sin(angle(sum((freqdiff ./ abs(freqdiff)),3)))) < sin(maxphaseerror)) = 0;
                    tmp(1:size(tmp,1)+1:end) = 0;
                    result.wplv(:,:,end+1) = tmp;
                    result.wplv_timestamp = [result.plv_timestamp timestamp{nstep}];
                    if isfield(cfg.CONNECT.PHASE,'storeangle') & cfg.CONNECT.PHASE.storeangle == true
                         tmp = abs(sin(angle(sum((freqdiff ./ abs(freqdiff)),3))));
                         tmp(1:size(tmp,1)+1:end) = 0;
                         result.wplv_angle(:,:,end+1) = tmp;
                    end
                case 'PLT'
                    phasediff = permute(phasediff,[1 3 2]);
                    plt = zeros(size(phasediff,1),size(phasediff,1));
                    for j = 1:size(phasediff,3)
                        phasediff2 = sign(sin(phasediff));
                        phasediff2(phasediff2 == 0) = 1;
                        phasediff2 = sum(diff(phasediff2,[],2)~=0,2);
                        plt(j,:) = 1 - exp(-(phasediff2+1).^-1*(size(phasediff,2)/header.samplingrate));
                        plt(j,j) = 0;
                    end
                    result.plt(:,:,end+1) = plt;
                    result.plt_timestamp = [result.plt_timestamp timestamp{nstep}];
                case 'PTE'
                    disp('      skip PTE, calculation for wavelet not possible');
            end
        end
        clearvars freqdiff phasediff
    end
    result.freqband = Freqs(Nfreq,3:4);
    result.PHASE.freqs = Freqs(Nfreq,3:4);
    result.PHASE.step = step;
    result.PHASE.firststep = lag+1;
    result.PHASE.laststep = (stepmax-1)*step + lag+1;
    if isfield(header,'channels')
        result.channels = header.channels;
    end
    
    if exist('Connect_file','var')
        % Store results
        if isfield(cfg,'patient')
            patient = cfg.patient; %#ok<NASGU>
        else
            patient = []; %#ok<NASGU>
        end
        save(fullfile(Connect_filepath,Connect_file),'result','patient');
    end
    Result.(Outputname) = result;
end