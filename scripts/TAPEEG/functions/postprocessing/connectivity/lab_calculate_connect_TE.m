% Calculate Synchronization likelyhood in specified frequency ranges
%
% [result,cfg] = lab_calculate_connect_TE(data,header,cfg)
%
% by F. Hatz 2012

function [Result,cfg,header] = lab_calculate_connect_TE(data,header,cfg)

if ~isfield(cfg,'Output_file')
    if exist('header','var') & isfield(header,'EEG_file')
        cfg.Output_file = header.EEG_file;
        cfg.Output_filepath = header.EEG_filepath;
    else
        cfg.Output_file = [];
    end
end
[~,~,~,cfg.Output_fileS] = lab_filename(cfg.Output_file);

if ~exist('lab_TE') %#ok<EXIST>
    Result = [];
    disp('    Abort: ''lab_TE'' missing')
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

if ~isfield(cfg.CONNECT,'TE') | ~isfield(cfg.CONNECT.TE,'el')
    cfg.CONNECT.TE.el = 5;
    cfg.CONNECT.TE.tau = 0.35;
end

Result = [];
IdxAll = [];
Epochs = [];
for Nfreq = 1:size(Freqs,1)
    headertmp = header;
    if max(Freqs(Nfreq,1:2)) == 0
        disp ('   Calculate TE unfiltered signal')
        Outputname = 'UNFILTERED';
    else
        disp (['   Calculate TE ' num2str(Freqs(Nfreq,3)) 'Hz to ' num2str(Freqs(Nfreq,4)) 'Hz'])
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
    
    % Calculate TE
    result.TE = [];
    [TE,TE_all,MI,MI_all,TE_delays] = lab_TE(datafilt,headertmp,cfg.CONNECT.TE,cfg);
    if isempty(TE)
        TE = NaN(size(datafilt,1),size(datafilt,1));
        TE_all = NaN(size(datafilt,1),size(datafilt,1));
        MI = NaN(size(datafilt,1),size(datafilt,1));
        MI_all = NaN(size(datafilt,1),size(datafilt,1));
        skippostprocessing = 1;
    else
        skippostprocessing = 0;
    end
    result.TE = TE;
    result.TE_all = TE_all;
    result.MI = MI;
    result.MI_all = MI_all;
    result.TE_delays = TE_delays(:)';
    
    result.freqband = Freqs(Nfreq,3:4);
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
        Connect_fileI = [Connect_file(1:end-4) '.block'];
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
        while exist(fullfile(Connect_filepath,Connect_fileI),'file')
            pause(3)
        end
        Fid = fopen(fullfile(Connect_filepath,Connect_fileI),'w');
        if exist(fullfile(Connect_filepath,Connect_file),'file')
            result2 = result;
            clearvars result
            try %#ok<TRYNC>
                load(fullfile(Connect_filepath,Connect_file));
            end
            if ~exist('result','var') | ~isfield(result,'TE')
                result.TE = [];
                result.TE_all = [];
                result.MI = [];
                result.MI_all = [];
                result.TE_delays = [];
            end
            result.TE = cat(3,result.TE,result2.TE);
            result.TE_all = cat(3,result.TE_all,result2.TE_all);
            result.MI = cat(3,result.MI,result2.MI);
            result.MI_all = cat(3,result.MI_all,result2.MI_all);
            result.TE_delays = cat(2,result.TE_delays,result2.TE_delays);
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
        fclose(Fid);
        delete(fullfile(Connect_filepath,Connect_fileI));
    end
    Result.(Outputname) = result;
    clearvars result
end