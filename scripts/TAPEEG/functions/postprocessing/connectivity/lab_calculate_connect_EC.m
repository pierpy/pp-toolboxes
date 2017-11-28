% Calculate envelope coherence
%
% AEC = average envelope correlation (correlation are calculated for every
%                                     step, average at the end)
% CAE = envelope average correlation (first average of step-data,
%                                     correlation-calculation at the end)
%
% [result,cfg] = lab_calculate_connect_EC(data,header,cfg)
%
% written by F. Hatz 2013

function [Result,cfg,header] = lab_calculate_connect_EC(data,header,cfg)

if isfield(cfg,'Output_file') & ~isempty(cfg.Output_file)
    [~,~,~,cfg.Output_fileS] = lab_filename(cfg.Output_file);
end

if ~exist('cfg','var') | ~isfield(cfg,'CONNECT') | ~isfield(cfg.CONNECT,'freqs')
    cfg.CONNECT.freqs = [1 30];
end
if size(cfg.CONNECT.freqs,2) == 4
    Freqs = cfg.CONNECT.freqs;
elseif size(cfg.CONNECT.freqs,2) == 2
    Freqs = [cfg.CONNECT.freqs cfg.CONNECT.freqs];
else
    Freqs = [0 0 0 0];
end

if ~isfield(cfg.CONNECT,'EC') | ~isfield(cfg.CONNECT.EC,'stepec')
    cfg.CONNECT.EC.stepec = 0.5;
end

Result = [];
IdxAll = [];
Epochs = [];
for Nfreq = 1:size(Freqs,1)
    headertmp = header;
    if max(Freqs(Nfreq,1:2)) == 0
        disp ('   Calculate envelope correlation - unfiltered signal')
        Outputname = 'UNFILTERED';
    else
        disp (['   Calculate envelope correlation ' num2str(Freqs(Nfreq,3)) 'Hz to ' num2str(Freqs(Nfreq,4)) 'Hz'])
        Outputname = ['F' num2str(Freqs(Nfreq,3)) 'F' num2str(Freqs(Nfreq,4))];
    end
    
    cfgfilt.lowpass = Freqs(Nfreq,2);
    cfgfilt.highpass = Freqs(Nfreq,1);
    if strcmp(cfg.CONNECT.filter,'freq-shift')
        disp('     set shifting for frequency filter')
        cfgfilt.filtermode = 'freq';
        cfgfilt.dofreqshift = true;
    else
        cfgfilt.filtermode = cfg.CONNECT.filter;
    end
    cfg.nodisp = 1;
    cfgfilt.nodisp = 1;
    
    % filter data
    if ~strcmp(cfg.CONNECT.filter,'input') 
        datafilt = lab_filter(data,headertmp,cfgfilt,'novrb');
    else
        datafilt = data;
    end
    
    % do hilbert if needed
    if isfield(cfg.CONNECT.EC,'EChilbert') & cfg.CONNECT.EC.EChilbert == true
        datafilt = hilbert(datafilt')';
        datafilt = abs(datafilt);
    else
        % take absolute values
        datafilt = abs(datafilt);
        % filter data
        cfgfilt = rmfield(cfgfilt,'highpass');
        datafilt = lab_filter(datafilt,headertmp,cfgfilt,'novrb');
    end
    
    if isempty(IdxAll)
        Idx = 1:size(datafilt,2);
        
        % exclude markers if necessary
        if isfield(cfg.CONNECT,'MARKER') & isfield(cfg.CONNECT.MARKER,'markerexclude') & ...
                iscell(cfg.CONNECT.MARKER.markerexclude) & ~isempty(cfg.CONNECT.MARKER.markerexclude) & ...
                isfield(headertmp,'events') & ~isempty(headertmp.events)
            disp('    Exclude periods with selected markers');
            [datafilt,headertmp,Idx2] = lab_exclude_markers(datafilt,headertmp,cfg.CONNECT.MARKER.markerexclude);
            if isempty(datafilt)
                disp('     Abort: calculation of connectivity not possible');
                return
            end
            Idx = Idx(Idx2);
            clearvars Idx2
        end
        
        % include markers if necessary
        if isfield(cfg.CONNECT,'MARKER') & isfield(cfg.CONNECT.MARKER,'markerinclude') & ...
                iscell(cfg.CONNECT.MARKER.markerinclude) & ~isempty(cfg.CONNECT.MARKER.markerinclude) & ...
                isfield(headertmp,'events') & ~isempty(headertmp.events)
            disp('    Restrict periods to selected markers');
            [datafilt,headertmp,Idx2] = lab_include_markers(datafilt,headertmp,cfg.CONNECT.MARKER.markerinclude);
            if isempty(datafilt)
                disp('     Abort: calculation of connectivity not possible');
                return
            end
            Idx = Idx(Idx2);
            clearvars Idx2
        end
        
        % Select TFs by minimal correlation (only for data with microstates)
        if isfield(headertmp,'CORR')
            headertmp.CORR = headertmp.CORR(:,Idx);
            if isfield(cfg.CONNECT,'MARKER') & isfield(cfg.CONNECT.MARKER,'mincorr') & cfg.CONNECT.MARKER.mincorr > 0
                disp(['     restict to minimal cluster-correlation of ' num2str(cfg.CONNECT.MARKER.mincorr)]);
                tmp = abs(headertmp.CORR) >= cfg.CONNECT.MARKER.mincorr;
                if ~isempty(tmp)
                    datafilt = datafilt(:,tmp);
                    headertmp.events = lab_reduce_events(headertmp.events,tmp);
                    headertmp.numtimeframes = size(datafilt,2);
                    clearvars tmp
                else
                    disp('     Abort: calculation of connectivity not possible');
                    return
                end
            end
        end
        IdxAll = Idx;
    else
        Idx = IdxAll;
        datafilt = datafilt(:,Idx);
        if isfield(headertmp,'CORR')
            headertmp.CORR = headertmp.CORR(:,Idx);
        end
        if isfield(headertmp,'events') & isfield(headertmp.events,'POS') & ~isempty(headertmp.events.POS)
            headertmp.events = lab_reduce_events(headertmp.events,Idx);
        end
    end
    
    if cfg.CONNECT.EC.stepec > 0
        step = floor(headertmp.samplingrate * cfg.CONNECT.EC.stepec);
    else
        step = floor(headertmp.samplingrate / min(Freqs(Nfreq,1:2)))*2;
    end
    if step < 100 % step should not drop below 100
        step = 100;
    end
    if step > size(datafilt,2)
        step = size(datafilt,2);
    end
    maxstep = floor(size(datafilt,2) / step);
    
    % find bad epochs
    if isfield(cfg.CONNECT,'EPOCH') & ~isempty(cfg.CONNECT.EPOCH) & maxstep > 1
        if isempty(Epochs) | size(Epochs.data,2) ~= step
            disp('     Find epochs with automated routines');
            cfg.CONNECT.EPOCH.length = step;
            cfg.CONNECT.EPOCH = lab_markerselect(cfg.CONNECT.EPOCH,cfg,headertmp);
            [Epochs,cfg.CONNECT.EPOCH] = lab_select_epochs(data(:,Idx),headertmp,cfg.CONNECT.EPOCH);
        end
        if ~isempty(Epochs)
            disp(['     ' num2str(size(Epochs.events.POS,2)) ' epochs found']);
            datatmp = zeros(size(datafilt,1),step,size(Epochs.events.POS,2));
            for Nepoch = 1:size(Epochs.events.POS,2)
                datatmp(:,:,Nepoch) = datafilt(:,Epochs.events.POS(Nepoch):(Epochs.events.POS(Nepoch)+step-1));
            end
            datafilt = datatmp;
            maxstep = size(datafilt,3);
            clearvars Nepoch datatmp
        else
            disp('      find epochs failed, process all data');
            datafilt = datafilt(:,1:step*maxstep);
            headertmp.numtimeframes = step*maxstep;
            datafilt = reshape(datafilt,size(datafilt,1),step,maxstep);
        end
    else
        datafilt = datafilt(:,1:step*maxstep);
        headertmp.numtimeframes = step*maxstep;
        datafilt = reshape(datafilt,size(datafilt,1),step,maxstep);
    end
    
    % do calculation
    resulttmp = zeros(size(datafilt,1),size(datafilt,1),maxstep);
    for i = 1:maxstep
        resulttmp(:,:,i) = corr(datafilt(:,:,i)');
    end
    result.AEC = abs(mean(resulttmp,3));
    result.CAE = abs(corr(mean(datafilt,3)'));
    result.EC.step = step;
    if isfield(headertmp,'channels')
        result.channels = headertmp.channels;
    end
    if isfield(headertmp,'timestamp')
        result.AEC_timestamp = cellstr([cfg.Output_fileS '_' num2str(headertmp.timestamp)]);
    elseif isfield(cfg,'Output_fileS')
        result.AEC_timestamp = cellstr(cfg.Output_fileS);
    else
        result.AEC_timestamp = cellstr('MeanEC');
    end
    result.CAE_timestamp = result.AEC_timestamp;
    
    % store results
    if isfield(cfg,'Output_filepath') & ~isempty(cfg.Output_filepath)
        Connect_filepath = cfg.Output_filepath;
        if isfield(cfg,'patient')
            Connect_file = [cfg.patient '_Conn_F' num2str(Freqs(Nfreq,3)) '_' num2str(Freqs(Nfreq,4)) '.mat'];
        else
            Connect_file = ['Conn_F' num2str(Freqs(Nfreq,3)) '_' num2str(Freqs(Nfreq,4)) '.mat'];
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
            result2 = result;
            clearvars result
            try %#ok<TRYNC>
                load(fullfile(Connect_filepath,Connect_file));
            end
            if ~exist('result','var') | ~isfield(result,'AEC')
                result.AEC = [];
                result.AEC_timestamp = [];
                result.CAE = [];
                result.CAE_timestamp = [];
            end
            result.AEC = cat(3,result.AEC,result2.AEC);
            result.CAE = cat(3,result.CAE,result2.CAE);
            result.AEC_timestamp = [result.AEC_timestamp result2.AEC_timestamp];
            result.CAE_timestamp = [result.CAE_timestamp result2.CAE_timestamp];
            result.EC.step = step;
        end
        patient = cfg.patient; %#ok<NASGU>
        if isfield(headertmp,'channels')
            result.channels = headertmp.channels;
        end
        if isfield(headertmp,'locs')
            result.locs = headertmp.locs;
        end
        result.freqband = Freqs(Nfreq,3:4);
        save(fullfile(Connect_filepath,Connect_file),'result','patient');
    end
    Result.(Outputname) = result;
    clearvars result
end