function [cfg,skipprocessing] = lab_set_extrachannels(cfg,header,data,skipdiag)

skipprocessing = 0;
clearExtra = 0;

if ~exist('skipdiag','var')
    skipdiag = 0;
end
if ~exist('data','var')
    data =[];
end
global THeader
if ~isempty(THeader)
    header = THeader;
elseif ~exist('header','var')
    header= [];
end
if ~exist('cfg','var')
    cfg = [];
end

if ~isfield(cfg,'EXTRA') | ~isfield(cfg.EXTRA,'numdatachans')
    if isfield(header,'numdatachannels')
        cfg.EXTRA.numdatachans = header.numdatachannels;
    elseif isfield(cfg,'EXTRA') & isfield(cfg.EXTRA,'ref_chan') & isnumeric(cfg.EXTRA.ref_chan) & length(cfg.MAIN.ref_chan) == 1
        cfg.EXTRA.numdatachans = cfg.EXTRA.ref_chan;
    else
        cfg.EXTRA.numdatachans = [];
        clearExtra = 1;
    end
end
if isfield(cfg.EXTRA,'numdatachans') & cfg.EXTRA.numdatachans == 116 & ~isfield(cfg.EXTRA,'reduceAAL')
    cfg.EXTRA.reduceAAL = false;
end
if ~isfield(cfg.EXTRA,'ref_chan')
    if isfield(header,'ref_chan') & ~isempty(header.ref_chan)
        if isnumeric(header.ref_chan)
            if isfield(header,'numchannels') & header.ref_chan <= header.numchannels & header.numchannels ~= 78
                cfg.EXTRA.ref_chan = header.ref_chan;
            else
                cfg.EXTRA.ref_chan = 'none';
            end
        else
            cfg.EXTRA.ref_chan=header.ref_chan;
        end
    elseif isfield(header,'numchannels') & header.numchannels == 78
        cfg.EXTRA.ref_chan = 'none';
    elseif isfield(header,'numdatachannels') == 257
        cfg.EXTRA.ref_chan = header.numdatachannels;
    else
        cfg.EXTRA.ref_chan = [];
        clearExtra = 1;
    end
end
if ~isfield(cfg.EXTRA,'auxmethod')
    if isfield(header,'numauxchannels') & ~isempty(header.numauxchannels) & header.numauxchannels > 0 & isfield(header,'numchannels')
        cfg.EXTRA.auxmethod = 'automatic';
    else
        cfg.EXTRA.auxmethod = 'none';
    end
    cfg.EXTRA.AUX = [];
end
if ~isfield(cfg.EXTRA,'ecg_ch')
    if isfield(header,'ecg_ch')
        cfg.EXTRA.ecg_ch = header.ecg_ch;
    else
        cfg.EXTRA.ecg_ch = [];
    end
end
if ~isfield(cfg.EXTRA,'computegfp')
    cfg.EXTRA.computegfp = false;
end

% show dialog
if isfield(cfg,'MAIN') & isfield(cfg.MAIN,'auto') & cfg.MAIN.auto == 1
    return;
end
if skipdiag == 1
    if clearExtra == 1
        cfg.EXTRA = [];
    end
    return
end

Prompt = cell(0,2);
Formats = [];

Prompt(end+1,:) = {'Number of data channels','numdatachans'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'integer';
Formats(end,1).limits = [0 999];
Formats(end,1).size = 40;
Formats(end,1).span = [1 4];

Prompt(end+1,:) = {'Reference channel (channelnumber(s) / none / string)',''};
Formats(end+1,1).type = 'text';
Formats(end,1).span = [1 4];

Prompt(end+1,:) = {'','ref_chan'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'text';
Formats(end,1).size = 80;
if isnumeric(cfg.EXTRA.ref_chan)
    cfg.EXTRA.ref_chan = num2str(cfg.EXTRA.ref_chan);
end
Formats(end,1).span = [1 4];

if isfield(cfg.EXTRA,'REFoverwrite') | isempty(data) | size(data,2) == 1 | ...
        (isnumeric(cfg.EXTRA.ref_chan) & max(cfg.EXTRA.ref_chan) <= size(data,1) & ...
        min(cfg.EXTRA.ref_chan) > 0 & mean(data(cfg.EXTRA.ref_chan,:)) ~= 0)
    Prompt(end+1,:) = {'Overwrite File-Reference','REFoverwrite'};
    Formats(end+1,1).type = 'check';
    Formats(end,1).span = [1 4];
    
    Prompt(end+1,:) = {'Subtract Reference-Channel from Data-Channels','REFcorrect'};
    Formats(end+1,1).type = 'check';
    Formats(end,1).span = [1 4];
else
    Formats(end+1,1).type = 'none';
end

Formats(end+1,1).type = 'none';
Formats(end,1).span = [1 4];

Prompt(end+1,:) = {'Auxillary channels','auxmethod'};
Formats(end+1,1).type = 'list';
Formats(end,1).style = 'popupmenu';
Formats(end,1).format = 'input';
Formats(end,1).items = {'automatic','none','channels'};
Formats(end,1).callback = {@set_auxmethod,'@ALL','@ALL'};
Formats(end,1).span = [1 2];

Prompt(end+1,:) = {'Channels','AUX'};
Formats(end+1,1).type = 'check';
Formats(end,1).callback = {@set_auxchans,'@ALL','@ALL'};

Formats(end+1,1).type = 'none';
Formats(end,1).span = [1 4];

Prompt(end+1,:) = {'Default ecg channel','ecg_ch'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'integer';
Formats(end,1).limits = [0 999];
Formats(end,1).size = 40;
Formats(end,1).span = [1 4];

Formats(end+1,1).type = 'none';
Formats(end,1).span = [1 4];

Prompt(end+1,:) = {'Compute GFP','computegfp'};
Formats(end+1,1).type = 'check';
Formats(end,1).size = [-1 -1];
Formats(end,1).span = [1 4];

if isfield(cfg.EXTRA,'removeNaN')
    Prompt(end+1,:) = {'Set invalid channels to zero','removeNAN'};
    Formats(end+1,1).type = 'check';
    Formats(end,1).span = [1 4];
end

if isfield(cfg.EXTRA,'reduceAAL')
    Prompt(end+1,:) = {'Reduce channels to 78','reduceAAL'};
    Formats(end+1,1).type = 'check';
    Formats(end,1).span = [1 4];
end

[cfg.EXTRA,Cancelled] = inputsdlg(Prompt,'Extra channels',Formats,cfg.EXTRA);
pause(0.1);
if Cancelled == 1
    cfg.EXTRA = [];
    skipprocessing = 1;
    return
elseif ~isnan(str2double(cfg.EXTRA.ref_chan))
    cfg.EXTRA.ref_chan = str2num(cfg.EXTRA.ref_chan); %#ok<ST2NM>
end
    
    function settings = set_auxmethod(settings)
        if strcmp(settings.auxmethod,'automatic') | strcmp(settings.auxmethod,'none')
            settings.AUX = [];
        else
            settings = set_auxchans(settings);
        end
    end
    
    function settings = set_auxchans(settings)
        settings.auxmethod = 'channels';
        if isfield(header,'locs') & ~isempty(header.locs)
            LOCS = header.locs;
        else
            LOCS = [];
        end
        if ~isempty(LOCS)
            if ~isempty(settings.AUX)
                LOCS.auxNr = settings.AUX.channels;
                LOCS.auxlabels = settings.AUX.labels;
                LOCS.aux = length(LOCS.auxNr);
            end
            settmp.LOCS = LOCS;
            settmp.Color = [1 1 1];
            settmp.ColorIdx = [1 0 0];
            settmp.indexed = [];
            settmp.Title = 'Set auxillary channels';
            [~,~,LOCS] = lab_plot_locs(settmp,1,0,0,0,1);
            if isfield(LOCS,'auxlabels') & isfield(LOCS,'auxNr') & ~isempty(LOCS.auxNr)
                settings.AUX.channels = LOCS.auxNr;
                settings.AUX.labels = LOCS.auxlabels;
            else
                settings.AUX = [];
                settings.auxmethod = 'automatic';
            end
        else
            if ~isempty(settings.AUX)
                AUX = cellstr(num2str(settings.AUX.channels'));
                AUX = cat(2,AUX,settings.AUX.labels');
            else
                AUX = {'',''};
            end
            if ~isempty(settings.numdatachans)
                Nchans = settings.numdatachans;
            else
                Nchans = 300;
            end
            AUX = lab_table_dialog(AUX,{'Channel','AuxName'},'Define AUX channels',1,{cellstr(num2str((1:Nchans)'))','char'});
            settings.AUX.channels = [];
            settings.AUX.labels = {};
            for i = 1:size(AUX,1)
                if ~isnan(str2double(AUX{i,1}))
                    settings.AUX.channels(1,end+1) = str2double(AUX{i,1});
                    settings.AUX.labels{1,end+1} = AUX{i,2};
                end
            end
            if isempty(settings.AUX.channels)
                settings.AUX = [];
            end
        end
        
    end
end