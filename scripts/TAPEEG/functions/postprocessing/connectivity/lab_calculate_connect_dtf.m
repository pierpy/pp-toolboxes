function [Result,cfg,header] = lab_calculate_connect_dtf(data,header,cfg)

if ~isfield(cfg,'Output_file')
    if exist('header','var') & isfield(header,'EEG_file')
        cfg.Output_file = header.EEG_file;
        cfg.Output_filepath = header.EEG_filepath;
    else
        cfg.Output_file = [];
    end
end
[~,~,~,cfg.Output_fileS] = lab_filename(cfg.Output_file);

if ~exist('cfg','var') | ~isfield(cfg,'CONNECT')
    cfg.CONNECT = [];
end
if ~isfield(cfg.CONNECT,'MARKER') | ~isfield(cfg.CONNECT.MARKER,'markerinclude')
    cfg.CONNECT = lab_get_MARKERS(cfg.CONNECT,data,header);
    if isempty(cfg.CONNECT.MARKER)
        Result = [];
        return
    end
end
if ~isfield(cfg.CONNECT,'DTF') | ~isfield(cfg.CONNECT.DTF,'lag')
    cfg.CONNECT = lab_get_DTF(cfg.CONNECT,data,header);
    if isempty(cfg.CONNECT.DTF)
        Result = [];
        return
    end
end

if ~exist('header','var')
    header = [];
end
if isfield(header,'numdatachannels')
    data = data(1:header.numdatachannels,:);
end

if ~exist('cfg','var') | ~isfield(cfg,'CONNECT') | ~isfield(cfg.CONNECT,'freqs') | isempty(cfg.CONNECT.freqs)
    if isfield(cfg.CONNECT.DTF,'lowfreq') & isfield(cfg.CONNECT.DTF,'highfreq')
        cfg.CONNECT.freqs = [cfg.CONNECT.DTF.lowfreq cfg.CONNECT.DTF.highfreq];
    else
        cfg.CONNECT.freqs = [1 30];
    end
end
if size(cfg.CONNECT.freqs,2) == 4
    Freqs = cfg.CONNECT.freqs;
elseif size(cfg.CONNECT.freqs,2) == 2
    Freqs = [cfg.CONNECT.freqs cfg.CONNECT.freqs];
else
    Freqs = [0 0 0 0];
end

% exclude markers if necessary
if isfield(cfg.CONNECT,'MARKER') & isfield(cfg.CONNECT.MARKER,'markerexclude') & ...
        iscell(cfg.CONNECT.MARKER.markerexclude) & ~isempty(cfg.CONNECT.MARKER.markerexclude) & ...
        isfield(header,'events') & ~isempty(header.events)
    disp('    Exclude periods with selected markers');
    [data,header] = lab_exclude_markers(data,header,cfg.CONNECT.MARKER.markerexclude);
    if isempty(data)
        disp('     Abort: calculation of connectivity not possible');
        return
    end
end

startpnt = [];
if isfield(cfg.CONNECT.MARKER,'markerinclude') & ~isempty(cfg.CONNECT.MARKER.markerinclude) & ...
        isfield(header,'events') & ~isempty(header.events)
    for i = 1:length(cfg.CONNECT.MARKER.markerinclude)
        if strcmp(cfg.CONNECT.MARKER.markerinclude{i},'Start')
            startpnt = cat(1,startpnt,1);
        else
            tmp = strcmp(header.events.TYP,cfg.CONNECT.MARKER.markerinclude{i});
            if ~isempty(tmp)
                startpnt = cat(1,startpnt,header.events.POS(1,tmp)');
            end
            clearvars tmp
        end
    end
end
if isempty(startpnt)
    startpnt = 1;
end
startpnt = startpnt((startpnt + cfg.CONNECT.DTF.length - 1) < size(data,2));
if isempty(startpnt)
    startpnt = 1;
    cfg.CONNECT.DTF.length = size(data,2);
end

Result = [];
for Nfreq = 1:size(Freqs,1)
    disp (['   Calculate DTF ' num2str(Freqs(Nfreq,3)) 'Hz to ' num2str(Freqs(Nfreq,4)) 'Hz'])
    Outputname = ['F' num2str(Freqs(Nfreq,3)) 'F' num2str(Freqs(Nfreq,4))];
    lowf = Freqs(Nfreq,2);
    highf = Freqs(Nfreq,1);
    
    result.DTF = zeros(size(data,1),size(data,1),length(startpnt));
    result.DTF_timestamp = {};
    for i = 1:length(startpnt)
        ts = data(:,startpnt(i):(startpnt(i)+cfg.CONNECT.DTF.length-1))';
        matrixtmp = DTF(ts,lowf,highf,cfg.CONNECT.DTF.order,header.samplingrate);
        if cfg.CONNECT.DTF.sigtest == true
            sig_dtfmatrix = DTFsigvalues(ts,lowf,highf,cfg.CONNECT.DTF.order, ...
                header.samplingrate,cfg.CONNECT.DTF.shufftimes,cfg.CONNECT.DTF.siglevel,[]);
            matrixtmp = DTFsigtest(matrixtmp,sig_dtfmatrix);
        end
        for j = 1:size(matrixtmp,1)
            matrixtmp(j,j,:) = 0;
        end
        scale = max(max(max(matrixtmp)));
        if scale == 0
            scale = 1;
        end
        result.DTF(:,:,i) = matrixtmp / scale;
        if length(startpnt) == 1
            if isfield(header,'timestamp')
                result.DTF_timestamp{1,end+1} = [cfg.Output_fileS '_' num2str(header.timestamp)];
            elseif isfield(cfg,'Output_fileS')
                result.DTF_timestamp{1,end+1} = cfg.Output_fileS;
            else
                result.DTF_timestamp{1,end+1} = 'DTF';
            end
        else
            if isfield(header,'timestamp')
                result.DTF_timestamp{1,end+1} = [cfg.Output_fileS '_' num2str(header.timestamp) '_' num2str(i)];
            elseif isfield(cfg,'Output_fileS')
                result.DTF_timestamp{1,end+1} = [cfg.Output_fileS '_' num2str(i)];
            else
                result.DTF_timestamp{1,end+1} = ['DTF' num2str(i)];
            end
        end
    end
    result.DTF_param = cfg.CONNECT.DTF;
    result.freqband = Freqs(Nfreq,3:4);

    % Store results
    if isfield(cfg,'Output_filepath') & ~isempty(cfg.Output_filepath)
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
            if ~exist('result','var') | ~isfield(result,'DTF')
                result.DTF = [];
                result.DTF_timestamp = [];
            end
            result.DTF = cat(3,result.DTF,result2.DTF);
            result.DTF_parm = result2.DTF_parm;
            result.DTF_timestamp = [result.DTF_timestamp result2.DTF_timestamp];
        end
        if isfield(header,'locs')
            result.locs = header.locs;
        end
        if isfield(header,'channels')
            result.channels = header.channels;
        end
        result.freqband = Freqs(Nfreq,3:4);
        patient = cfg.patient; %#ok<NASGU>
        save(fullfile(Connect_filepath,Connect_file),'result','patient');
    end
    Result.(Outputname) = result;
    clearvars result
end

end