% Calculate Synchronization likelyhood in specified frequency ranges
%
% [result,cfg] = lab_calculate_connect_SL(data,header,cfg)
%
% by F. Hatz 2012

function [Result,cfg,header] = lab_calculate_connect_SL(data,header,cfg)

if ~isfield(cfg,'Output_file')
    if exist('header','var') & isfield(header,'EEG_file')
        cfg.Output_file = header.EEG_file;
        cfg.Output_filepath = header.EEG_filepath;
    else
        cfg.Output_file = [];
    end
end
[~,~,~,cfg.Output_fileS] = lab_filename(cfg.Output_file);

if ~exist('lab_SL') %#ok<EXIST>
    Result = [];
    disp('    Abort: ''lab_SL'' missing')
    return
end

if ~exist('cfg','var') | ~isfield(cfg,'CONNECT') | ~isfield(cfg.CONNECT,'freqs')
    cfg.CONNECT.freqs = [1 30];
end
if size(cfg.CONNECT.freqs,2) == 4
    Freqs = cfg.CONNECT.freqs;
elseif size(cfg.CONNECT.freqs,2) == 2
    Freqs = [cfg.CONNECT.freqs cfg.CONNECT.freqs];
else
    Freqs = 0;
end

if ~isfield(cfg.CONNECT,'SL') | ~isfield(cfg.CONNECT.SL,'step')
    cfg.CONNECT.SL.step = 0;
end
if ~isfield(cfg.CONNECT.SL,'pref')
    if ~isfield(cfg,'MAIN') | ~isfield(cfg.MAIN,'auto') | cfg.MAIN.auto == 0
        disp ('   Calculate SL: Ask for pref')
        prompt={'pref'};
        name='pref';
        numlines(1,1) = ones(1,1);
        numlines(1,2) = ones(1,1)*20;
        if isfield(cfg.CONNECT.SL,'pref')
            defaultanswer={num2str(cfg.CONNECT.SL.pref)};
        else
            defaultanswer={'0.05'};
        end
        answer=inputdlg(prompt,name,numlines,defaultanswer);
        pause(0.1)
        if size(answer,1) > 0
            cfg.CONNECT.SL.pref = str2num(answer{1,1}); %#ok<ST2NM>
        else
            cfg.CONNECT.SL.pref = 0.05;
            disp('   SL-pref set to 0.05')
        end
        clearvars defaultanswer answer name numlines prompt tmp
    else
        cfg.CONNECT.SL.pref = 0.05;
        disp('    SL-pref set to 0.05')
    end
end

if ~isfield(cfg.CONNECT.SL,'storehits')
    cfg.CONNECT.SL.storehits = false;
end

Result = [];
IdxAll = [];
Epochs = [];
for Nfreq = 1:size(Freqs,1)
    headertmp = header;
    if max(Freqs(Nfreq,1:2)) == 0
        disp ('   Calculate SL unfiltered signal')
        Outputname = 'UNFILTERED';
    else
        disp (['   Calculate SL ' num2str(Freqs(Nfreq,3)) 'Hz to ' num2str(Freqs(Nfreq,4)) 'Hz'])
        Outputname = ['F' num2str(Freqs(Nfreq,3)) 'F' num2str(Freqs(Nfreq,4))];
    end
    
    % Filter data
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
    if ~strcmp(cfgfilt.filtermode,'input')
        datafilt = lab_filter(data,headertmp,cfgfilt,'novrb');
    else
        datafilt = data;
    end
    
    if isempty(IdxAll)
        Idx = 1:size(datafilt,2);
        
        % exclude markers if necessary
        if isfield(cfg.CONNECT,'MARKER') & isfield(cfg.CONNECT.MARKER,'markerexclude') & ...
                iscell(cfg.CONNECT.MARKER.markerexclude) & ~isempty(cfg.CONNECT.MARKER.markerexclude) & ...
                isfield(headertmp,'events') & ~isempty(headertmp.events)
            disp('    Exclude periods with selected markers');
            [datafilt,headertmp,Idx2] = lab_exclude_markers(datafilt,headertmp,cfg.CONNECT.MARKER.markerexclude);
            Idx = Idx(1,Idx2);
            clearvars Idx2
        end
        
        % include markers if necessary
        if isfield(cfg.CONNECT,'MARKER') & isfield(cfg.CONNECT.MARKER,'markerinclude') & ...
                iscell(cfg.CONNECT.MARKER.markerinclude) & ~isempty(cfg.CONNECT.MARKER.markerinclude) & ...
                isfield(headertmp,'events') & ~isempty(headertmp.events)
            disp('    Restrict periods to selected markers');
            [datafilt,headertmp,Idx2] = lab_include_markers(datafilt,headertmp,cfg.CONNECT.MARKER.markerinclude);
            Idx = Idx(1,Idx2);
            clearvars Idx2
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
    
    % set parameters
    speed = floor(headertmp.samplingrate / max(Freqs(Nfreq,1:2)));
    if speed < 60 % step should not drop below 60
        speed = 60;
    end
    if isempty(speed)
        speed = 100;
    end
    lag = round(headertmp.samplingrate / (3*Freqs(Nfreq,2)));
    M = round(3 * (Freqs(Nfreq,2) / Freqs(Nfreq,1)))+1;
    W1 = 2 * lag * (M - 1);
    pref = cfg.CONNECT.SL.pref;
    if isfield(cfg.CONNECT.SL,'step') & cfg.CONNECT.SL.step > 1
        W2 = 10/pref + W1 - 1;
    else
        W2 = size(datafilt,2);
    end
    if pref*(W2-W1+1) < 10
        pref = 10/(W2-W1+1);
    end
    
    % calculate step
    if cfg.CONNECT.SL.step > 0 & strcmp(cfg.CONNECT.SL.stepunit,'TFs') & cfg.CONNECT.SL.step <= size(datafilt,2);
        step = cfg.CONNECT.SL.step;
    elseif cfg.CONNECT.SL.step > 0 & strcmp(cfg.CONNECT.SL.stepunit,'seconds') & ...
            isfield(headertmp,'samplingrate') & cfg.CONNECT.SL.step*headertmp.samplingrate <= size(datafilt,2);
        step = cfg.CONNECT.SL.step * headertmp.samplingrate;
    elseif cfg.CONNECT.SL.step > 0 & strcmp(cfg.CONNECT.SL.stepunit,'phases of max frequency')
        step = floor(headertmp.samplingrate * (cfg.CONNECT.SL.step / max(Freqs(Nfreq,1:2))));
    else
        step = size(datafilt,2);
    end
    if step < 60 % step should not drop below 60
        step = 60;
    end
    if step > size(datafilt,2);
        step = size(datafilt,2);
    end
    if W2 > size(datafilt,2);
        W2 = size(datafilt,2);
    end
    stepmax = floor((size(datafilt,2) - W2 + step) / step);
    if stepmax < 1
        stepmax = 1;
        step = size(datafilt,2);
    end
    if stepmax > 1
        cfg.CONNECT.SL.storesingle = true;
    end
    
    % Find bad epochs
    if isfield(cfg.CONNECT,'EPOCH') & ~isempty(cfg.CONNECT.EPOCH) & stepmax > 1
        if isempty(Epochs) | size(Epochs.data,2) ~= step
            disp('     Find epochs with automated routines');
            cfg.CONNECT.EPOCH.length = step;
            cfg.CONNECT.EPOCH = lab_markerselect(cfg.CONNECT.EPOCH,cfg,headertmp);
            [Epochs,cfg.CONNECT.EPOCH] = lab_select_epochs(data(:,Idx),headertmp,cfg.CONNECT.EPOCH,W2);
        end
        if ~isempty(Epochs)
            disp(['     ' num2str(size(Epochs.events.POS,2)) ' epochs found']);
            datatmp = zeros(size(datafilt,1),step,size(Epochs.events.POS,2));
            for Nepoch = 1:size(Epochs.events.POS,2)
                datatmp(:,:,Nepoch) = datafilt(:,Epochs.events.POS(Nepoch):(Epochs.events.POS(Nepoch)+step-1));
            end
            datafilt = datatmp;
            stepmax = size(datafilt,3);
            clearvars Nepoch datatmp
        else
            disp('      find epochs failed, process all data');
            datatmp = zeros(size(datafilt,1),W2,stepmax);
            for nstep = 1:stepmax
                datatmp(:,:,nstep) = datafilt(:,(nstep-1)*step + 1:(nstep-1)*step + W2);
            end
            datafilt = datatmp;
            clearvars datatmp
        end
    else
        datafilt = datafilt(:,1:(stepmax-1)*step+W2);
        datatmp = zeros(size(datafilt,1),W2,stepmax);
        for nstep = 1:stepmax
            datatmp(:,:,nstep) = datafilt(:,(nstep-1)*step + 1:(nstep-1)*step + W2);
        end
        datafilt = datatmp;
        clearvars datatmp
    end
    
    % Calculate SL
    result.SL = [];
    result.SL_hit = [];
    result.SL_timestamp = {};
    if stepmax > 1
        fprintf('   SL')
    end
    for nstep = 1:stepmax
        if stepmax > 1
            fprintf('.')
        end
        datatmp = datafilt(:,:,nstep);
        [SL,SL_hit] = lab_SL(datatmp',lag,M,W1,W2,pref,speed,'s',cfg.CONNECT.SL.storehits);
        if isempty(SL)
            SL = NaN(size(datafilt,1),size(datafilt,1));
            SL_hit = SL;
            skippostprocessing = 1;
        else
            skippostprocessing = 0;
        end
        if skippostprocessing == 0
            if isfield(cfg.CONNECT.SL,'storesingle') & cfg.CONNECT.SL.storesingle == true
                result.SL = cat(3,result.SL,mean(SL,3));
                result.SL_hit = cat(3,result.SL_hit,SL_hit);
                if isfield(headertmp,'timestamp')
                    result.SL_timestamp{1,end+1} = [cfg.Output_fileS '_' num2str(headertmp.timestamp)];
                elseif isfield(cfg,'Output_fileS')
                    result.SL_timestamp{1,end+1} = cfg.Output_fileS;
                else
                    result.SL_timestamp{1,end+1} = 'MeanSL';
                end
            else
                result.SL = cat(3,result.SL,SL);
                result.SL_hit = cat(3,result.SL_hit,SL_hit);
                for i = 1:size(result.SL,3)
                    if isfield(headertmp,'timestamp')
                        result.SL_timestamp{1,end+1} = [cfg.Output_fileS '_' num2str(headertmp.timestamp) '_' num2str(i)];
                    elseif isfield(cfg,'Output_fileS')
                        result.SL_timestamp{1,end+1} = [cfg.Output_fileS '_' num2str(i)];
                    else
                        result.SL_timestamp{1,end+1} = ['SL' num2str(i)];
                    end
                end
            end
        end
    end
    if stepmax > 1
        disp(':')
    end
    result.freqband = Freqs(Nfreq,3:4);
    result.SL_parm.speed = speed;
    result.SL_parm.lag = lag;
    result.SL_parm.M = M;
    result.SL_parm.W1 = W1;
    result.SL_parm.W2 = W2;
    result.SL_parm.pref = pref;
    if isfield(headertmp,'channels')
        result.channels = headertmp.channels;
    end
    
    % Store results
    if isfield(cfg,'Output_filepath') & ~isempty(cfg.Output_filepath) & ...
            skippostprocessing == 0
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
            elseif min(~strcmp(cfg.listold,fullfile(Connect_filepath,Connect_file)))
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
            if ~exist('result','var') | ~isfield(result,'SL')
                result.SL = [];
                result.SL_hit = [];
                result.SL_timestamp = [];
            end
            result.SL = cat(3,result.SL,result2.SL);
            if ~strcmp(mexext,'mexw32') &  cfg.CONNECT.SL.storehits == true
                result.SL_hit = cat(3,result.SL_hit,result2.SL_hit);
            end
            result.SL_parm =result2.SL_parm;
            result.SL_timestamp = [result.SL_timestamp result2.SL_timestamp];
        end
        if isfield(headertmp,'locs')
            result.locs = headertmp.locs;
        end
        if isfield(headertmp,'channels')
            result.channels = headertmp.channels;
        end
        result.freqband = Freqs(Nfreq,3:4);
        patient = cfg.patient; %#ok<NASGU>
        save(fullfile(Connect_filepath,Connect_file),'result','patient');
    end
    Result.(Outputname) = result;
    clearvars result
end