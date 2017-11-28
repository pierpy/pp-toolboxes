% Calculate Synchronization likelyhood in specified frequency ranges
% with vcorrection for volume conduction (using linear regression)
%
% [result,cfg] = lab_calculate_connect_SLc(data,header,cfg)
%
% by F. Hatz 2012

function [Result,cfg,header] = lab_calculate_connect_SLc(data,header,cfg)

if ~isfield(cfg,'Output_file')
    if exist('header','var') & isfield(header,'EEG_file')
        cfg.Output_file = header.EEG_file;
        cfg.Output_filepath = header.EEG_filepath;
    else
        cfg.Output_file = [];
    end
end
[~,~,~,cfg.Output_fileS] = lab_filename(cfg.Output_file);

if ~exist('lab_SLcorrect') %#ok<EXIST>
    Result = [];
    disp('    Abort: ''lab_SLcorrect'' missing')
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
    Freqs = [0 0 0 0];
end

if ~isfield(cfg.CONNECT,'SLc') | ~isfield(cfg.CONNECT.SLc,'step')
    cfg.CONNECT.SLc.step = 0;
end
if ~isfield(cfg.CONNECT.SLc,'pref')
    if ~isfield(cfg,'MAIN') | ~isfield(cfg.MAIN,'auto') | cfg.MAIN.auto == 0
        disp ('   Calculate SLc: Ask for pref')
        prompt={'pref'};
        name='pref';
        numlines(1,1) = ones(1,1);
        numlines(1,2) = ones(1,1)*20;
        if isfield(cfg.CONNECT.SLc,'pref')
            defaultanswer={num2str(cfg.CONNECT.SLc.pref)};
        else
            defaultanswer={'0.01'};
        end
        answer=inputdlg(prompt,name,numlines,defaultanswer);
        pause(0.1)
        if size(answer,1) > 0
            cfg.CONNECT.SLc.pref=str2num(answer{1,1});
        else
            cfg.CONNECT.SLc.pref = 0.01;
            disp('   SLc-pref set to 0.01')
        end
        clearvars defaultanswer answer name numlines prompt tmp
    else
        cfg.CONNECT.SLc.pref = 0.01;
        disp('    SLc-pref set to 0.01')
    end
end

if ~isfield(cfg.CONNECT.SLc,'storehits')
    cfg.CONNECT.SLc.storehits = false;
end
if ~isfield(cfg.CONNECT.SLc,'Mcorrect')
    cfg.CONNECT.SLc.Mcorrect = 1;
end
Result = [];
IdxAll = [];
Epochs = [];
for Nfreq = 1:size(Freqs,1)
    headertmp = header;
    if max(Freqs(Nfreq,:)) == 0
        disp ('   Calculate SLc unfiltered signal')
        Outputname = 'UNFILTERED';
    else
        disp (['   Calculate SLc ' num2str(Freqs(Nfreq,3)) 'Hz to ' num2str(Freqs(Nfreq,4)) 'Hz'])
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
    pref = cfg.CONNECT.SLc.pref;
    if isfield(cfg.CONNECT.SLc,'step') & cfg.CONNECT.SLc.step > 1
        W2 = 10/pref + W1 - 1;
    else
        W2 = size(datafilt,2);
    end
    if pref * (W2-W1+1) < 10
        pref = 10/(W2-W1+1);
    end
    
    % calculate step
    if cfg.CONNECT.SLc.step > 0 & strcmp(cfg.CONNECT.SLc.stepunit,'TFs') & cfg.CONNECT.SLc.step <= size(datafilt,2);
        step = cfg.CONNECT.SLc.step;
    elseif cfg.CONNECT.SLc.step > 0 & strcmp(cfg.CONNECT.SLc.stepunit,'seconds') & ...
            isfield(headertmp,'samplingrate') & cfg.CONNECT.SLc.step*headertmp.samplingrate <= size(datafilt,2);
        step = cfg.CONNECT.SLc.step * headertmp.samplingrate;
    elseif cfg.CONNECT.SLc.step > 0 & strcmp(cfg.CONNECT.SLc.stepunit,'phases of max frequency')
        step = floor(headertmp.samplingrate * (cfg.CONNECT.SLc.step / max(Freqs(Nfreq,1:2))));
    else
        step = size(datafilt,2);
    end
    if step < 60 % step should not drop below 60
        step = 60;
    end
    if step > size(datafilt,2);
        step = size(datafilt,2);
    end
    stepmax = floor((size(datafilt,2) - W2 + step) / step);
    if stepmax < 1
        stepmax = 1;
        step = size(datafilt,2);
    end
    if stepmax > 1
        cfg.CONNECT.SLc.storesingle = true;
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
    
    result.SLc = [];
    result.SLc_hit = [];
    result.SLc_timestamp = {};
    for nstep = 1:stepmax
        if stepmax > 1
            disp(['      step ' num2str(nstep)])
        end
        datatmp = datafilt(:,:,nstep);
        [SLc,SLc_hit] = lab_SLcorrect(datatmp',lag,M,W1,W2,pref,speed,'s',cfg.CONNECT.SLc.storehits,cfg.CONNECT.SLc.Mcorrect);
        if isempty(SLc)
            SLc = NaN(size(datafilt,1),size(datafilt,1));
            SLc_hit = SLc;
            skippostprocessing = 1;
        else
            skippostprocessing = 0;
        end
        if skippostprocessing == 0
            if isfield(cfg.CONNECT.SLc,'storesingle') & cfg.CONNECT.SLc.storesingle == true
                result.SLc = cat(3,result.SLc,mean(SLc,3));
                result.SLc_hit = cat(3,result.SLc_hit,SLc_hit);
                if isfield(headertmp,'timestamp')
                    result.SLc_timestamp{1,end+1} = [cfg.Output_fileS '_' num2str(headertmp.timestamp)];
                elseif isfield(cfg,'Output_fileS')
                    result.SLc_timestamp{1,end+1} = cfg.Output_fileS;
                else
                    result.SLc_timestamp{1,end+1} = 'MeanSLc';
                end
            else
                result.SLc = cat(3,result.SLc,SLc);
                result.SLc_hit = cat(3,result.SLc_hit,SLc_hit);
                for i = 1:size(result.SLc,3)
                    if isfield(headertmp,'timestamp')
                        result.SLc_timestamp{1,end+1} = [cfg.Output_fileS '_' num2str(headertmp.timestamp) '_' num2str(i)];
                    elseif isfield(cfg,'Output_fileS')
                        result.SLc_timestamp{1,end+1} = [cfg.Output_fileS '_' num2str(i)];
                    else
                        result.SLc_timestamp{1,end+1} = ['SLc' num2str(i)];
                    end
                end
            end
        end
    end
    result.freqband = Freqs(Nfreq,3:4);
    result.SLc_parm.speed = speed;
    result.SLc_parm.lag = lag;
    result.SLc_parm.M = M;
    result.SLc_parm.W1 = W1;
    result.SLc_parm.W2 = W2;
    result.SLc_parm.pref = pref;
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
            if ~exist('result','var') | ~isfield(result,'SLc')
                result.SLc = [];
                result.SLc_hit = [];
                result.SLc_timestamp = [];
            end
            result.SLc = cat(3,result.SLc,result2.SLc);
            if ~strcmp(mexext,'mexw32') &  cfg.CONNECT.SLc.storehits == true
                result.SLc_hit = cat(3,result.SLc_hit,result2.SLc_hit);
            end
            result.SLc_parm =result2.SLc_parm;
            result.SLc_timestamp = [result.SLc_timestamp result2.SLc_timestamp];
        end
        if isfield(headertmp,'locs')
            result.locs = headertmp.locs;
        end
        if isfield(headertmp,'channels')
            result.channels = headertmp.channels;
        end
        patient = cfg.patient; %#ok<NASGU>
        save(fullfile(Connect_filepath,Connect_file),'result','patient');
    end
    Result.(Outputname) = result;
    clearvars result
end