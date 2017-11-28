function [cfg,skipprocessing] = lab_set_filtering(cfg,header)

skipprocessing = 0;

global THeader
if ~isempty(THeader)
    header = THeader;
elseif ~exist('header','var')
    header = [];
end

if ~exist('cfg','var')
    cfg = [];
end

if ~isfield(header,'numchannels')
    header.numchannels = 300;
end
if ~isfield(header,'numdatachannels')
    header.numdatachannels = header.numchannels;
end
if isfield(header,'channels')
    Channels = cellstr(header.channels);
else
    Channels = cellstr((num2str((1:header.numchannels)')));
end

if ~exist('simple','var')
    disp ('   Ask for filter-settings')
end

if ~isfield(cfg,'FILT') | isempty(cfg.FILT)
    cfg.FILT.filtermode = 'firls';
    cfg.FILT.filtorder = 0;
    cfg.FILT.wavsmoothing = [];
    cfg.FILT.filterdefs = {};
    cfg.FILT.filterdefs{1,1} = Channels{1};
    cfg.FILT.filterdefs{1,2} = Channels{header.numdatachannels};
    cfg.FILT.filterdefs{1,3} = 50;
    cfg.FILT.filterdefs{1,4} = 1;
    cfg.FILT.filterdefs{1,5} = 70;
    if header.numdatachannels < header.numchannels
        cfg.FILT.filterdefs{2,1} = Channels{header.numdatachannels+1};
        cfg.FILT.filterdefs{2,2} = Channels{header.numchannels};
        cfg.FILT.filterdefs{2,3} = 50;
        cfg.FILT.filterdefs{2,4} = 1;
        cfg.FILT.filterdefs{2,5} = 70;
    end
else
    if isfield(cfg.FILT,'Filter')
        cfg.FILT.filterdefs = {};
        try %#ok<TRYNC>
            for i = 1:length(cfg.FILT.Filter)
                if cfg.FILT.Filter(i).firstchan >= 1 & cfg.FILT.Filter(i).firstchan <= header.numchannels & ...
                        cfg.FILT.Filter(i).lastchan >= 1 & cfg.FILT.Filter(i).lastchan <= header.numchannels
                    cfg.FILT.filterdefs{end+1,1} = Channels{cfg.FILT.Filter(i).firstchan};
                    cfg.FILT.filterdefs{end,2} = Channels{cfg.FILT.Filter(i).lastchan};
                    cfg.FILT.filterdefs{end,3} = cfg.FILT.Filter(i).notch;
                    cfg.FILT.filterdefs{end,4} = cfg.FILT.Filter(i).highpass;
                    cfg.FILT.filterdefs{end,5} = cfg.FILT.Filter(i).lowpass;
                else
                    disp('     skip original filter-setting, invalid first or last channel')
                end
            end
        end
        if isempty(cfg.FILT.filterdefs)
            cfg.FILT.filtermode = 'firls';
            cfg.FILT.filtorder = 0;
            cfg.FILT.wavsmoothing = [];
            cfg.FILT.filterdefs = {};
            cfg.FILT.filterdefs{1,1} = Channels{1};
            cfg.FILT.filterdefs{1,2} = Channels{header.numdatachannels};
            cfg.FILT.filterdefs{1,3} = 50;
            cfg.FILT.filterdefs{1,4} = 1;
            cfg.FILT.filterdefs{1,5} = 70;
            if header.numdatachannels < header.numchannels
                cfg.FILT.filterdefs{2,1} = Channels{header.numdatachannels+1};
                cfg.FILT.filterdefs{2,2} = Channels{header.numchannels};
                cfg.FILT.filterdefs{2,3} = 50;
                cfg.FILT.filterdefs{2,4} = 1;
                cfg.FILT.filterdefs{2,5} = 70;
            end
        end
    elseif isfield(cfg.FILT,'highpass') % compatibility for old settings-file
        cfg.FILT.filterdefs{1,1} = Channels{1};
        cfg.FILT.filterdefs{1,2} = Channels{header.numdatachannels};
        cfg.FILT.filterdefs{1,3} = cfg.FILT.notch;
        cfg.FILT.filterdefs{1,4} = cfg.FILT.highpass;
        cfg.FILT.filterdefs{1,5} = cfg.FILT.lowpass;
        if isfield(cfg.FILT,'filterauxchannels ') & cfg.FILT.filterauxchannels == true & header.numdatachannels < header.numchannels
            cfg.FILT.filterdefs{2,1} = Channels{header.numdatachannels+1};
            cfg.FILT.filterdefs{2,2} = Channels{header.numchannels};
            cfg.FILT.filterdefs{2,3} = cfg.FILT.notch;
            cfg.FILT.filterdefs{2,4} = cfg.FILT.highpass;
            cfg.FILT.filterdefs{2,5} = cfg.FILT.lowpass;
        end
    end
end

if strcmp(cfg.FILT.filtermode,'freq') & isfield(cfg.FILT,'dofreqshift') & cfg.FILT.dofreqshift == 1
    cfg.FILT.filtermode = 'freq-shift';
end
if isfield(cfg.FILT,'off')
    cfg.FILT = rmfield(cfg.FILT,'off');
end

Prompt = cell(0,2);
Formats = [];

Prompt(end+1,:) = {'Filter','filtermode'};
Formats(end+1,1).type = 'list';
Formats(end,1).style = 'popupmenu';
Formats(end,1).format = 'input';
Formats(end,1).items = {'firls','butter','freq','wavelet','cheby','freq-shift'};
Formats(end,1).callback = {@lab_get_filtermode,'@ALL','@ALL',1};

Prompt(end+1,:) = {'Order','filtorder'};
Formats(end+1,1).type = 'check';
Formats(end,1).callback = {@lab_get_filtermode,'@ALL','@ALL'};

Prompt(end+1,:) = {'WaveletSmoothing','wavsmoothing'};
Formats(end+1,1).type = 'check';
Formats(end,1).callback = {@lab_get_filtermode,'@ALL','@ALL'};

Formats(end+1,1).type = 'none';
Formats(end+1,1).type = 'none';

Prompt(end+1,:) = {'','filterdefs'};
Formats(end+1,1).type = 'table';
Formats(end,1).format = 'table';
Formats(end,1).items = {{'First channel','Last channel','Notch','Highpass','Lowpass'},{},{},{Channels(:)',Channels(:)','numeric','numeric','numeric'}};
Formats(end,1).size = [430 130];
Formats(end,1).span = [1 5];

Prompt(end+1,:) = {'Add filter','Button'};
Formats(end+1,1).type = 'button';
Formats(end,1).style = 'pushbutton';
Formats(end,1).size = [80 25];
Formats(end,1).callback = {@add_filter,'filterdefs','filterdefs'};

Prompt(end+1,:) = {'Delete filter','Button'};
Formats(end+1,1).type = 'button';
Formats(end,1).style = 'pushbutton';
Formats(end,1).size = [80 25];
Formats(end,1).callback = {@delete_filter,'filterdefs','filterdefs'};

[cfg.FILT,Cancelled] = inputsdlg(Prompt,'Filter',Formats,cfg.FILT);
if isempty(cfg.FILT) | Cancelled == 1
    cfg.FILT = [];
    skipprocessing = 1;
else
    pause(0.2);
    if strcmp(cfg.FILT.filtermode,'freq-shift')
        cfg.FILT.filtermode = 'freq';
        cfg.FILT.dofreqshift = 1;
    else
        cfg.FILT.dofreqshift = 0;
    end
    cfg.FILT.Filter = [];
    FlagChans = false(1,length(Channels));
    for i = 1:size(cfg.FILT.filterdefs,1)
        if isempty(cfg.FILT.filterdefs{i,1}) | isempty(find(strcmp(Channels,cfg.FILT.filterdefs{i,1}),1,'first'))
            break
        end
        Firstchan = find(strcmp(Channels,cfg.FILT.filterdefs{i,1}),1,'first');
        Lastchan = find(strcmp(Channels,cfg.FILT.filterdefs{i,2}),1,'first');
        tmp = find(FlagChans(Firstchan:Lastchan)==0,1,'first');
        if isempty(tmp)
            break
        end
        tmp2 = find(FlagChans(Firstchan+tmp:Lastchan)==1,1,'first');
        if isempty(tmp2)
            tmp2 = Lastchan - Firstchan + 1;
        else
            tmp2 = tmp2 + tmp - 1;
        end
        if tmp2 <= tmp
            break
        end
        Lastchan = Firstchan + tmp2 - 1;
        Firstchan = Firstchan + tmp - 1;
        cfg.FILT.Filter(end+1).firstchan = Firstchan;
        cfg.FILT.Filter(end).lastchan = Lastchan;
        if ~isnan(cfg.FILT.filterdefs{i,3})
            cfg.FILT.Filter(end).notch = cfg.FILT.filterdefs{i,3};
        else
            cfg.FILT.Filter(end).notch = [];
        end
        if ~isnan(cfg.FILT.filterdefs{i,4})
            cfg.FILT.Filter(end).highpass = cfg.FILT.filterdefs{i,4};
        else
            cfg.FILT.Filter(end).highpass = [];
        end
        if ~isnan(cfg.FILT.filterdefs{i,5})
            cfg.FILT.Filter(end).lowpass = cfg.FILT.filterdefs{i,5};
        else
            cfg.FILT.Filter(end).lowpass = [];
        end
        FlagChans(Firstchan:Lastchan) = true;
    end
    cfg.FILT = rmfield(cfg.FILT,'filterdefs');
end

    function filterdefs = add_filter(filterdefs)
        filterdefs{end+1,1} = '';
        filterdefs{end,2} = '';
    end
    
    function filterdefs = delete_filter(filterdefs)
        if size(filterdefs,1) > 1
            filterdefs = filterdefs(1:end-1,:);
        else
            filterdefs = {};
            filterdefs{1,1} = Channels{1};
            filterdefs{1,2} = Channels{header.numdatachannels};
            filterdefs{1,3} = 50;
            filterdefs{1,4} = 1;
            filterdefs{1,5} = 70;
        end
    end
    
end