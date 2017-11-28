% Load eeg/meg-file for manual processing
%
% Written by F. Hatz 2013 Neurology Basel

function data = lab_get_data(data)

if ~exist('data','var')
    data = [];
end
header = lab_create_header(data,true);
cfg = [];   

settings.Filename = '';
settings.data = header2string(header);

Prompt = cell(0,2);
Formats = [];

Prompt(end+1,:) = {'Filename',''};
Formats(end+1,1).type = 'text';
Formats(end,1).span = [1 2];

Prompt(end+1,:) = {'','Filename'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'file';
Formats(end,1).items = {'*.*','EEG/MEG-File'};
Formats(end,1).limits = [0 1]; % single file get
Formats(end,1).callback =  {@read_eeg,'@ALL','@ALL'};
Formats(end,1).size = 350;
Formats(end,1).span = [1 2];

Prompt(end+1,:) = {'','data'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'text';
Formats(end,1).limits = [0 9];
Formats(end,1).span = [5,1];
Formats(end,1).size = [350 100];
Formats(end,1).enable = 'inactive';
Formats(end,1).span = [1 2];

Prompt(end+1,:) = {'Edit',''};
Formats(end+1,1).type = 'button';
Formats(end,1).style = 'pushbutton';
Formats(end,1).size = [100 25];
Formats(end,1).callback = {@do_edit,'@ALL','@ALL'};

[~,Cancelled] = inputsdlg(Prompt,'EEG/MEG',Formats,settings);
if Cancelled == 1
    data = [];
elseif isfield(header,'numdatachannels')
    data = data(1:header.numdatachannels,:);
end

% nested functions
    function settings = read_eeg(settings)
        if exist(settings.Filename,'file')
            Filename = settings.Filename;
        else
            Filename = '';
        end
        [data,header,cfg] = lab_read_data(Filename,cfg);
        if ~isempty(data) & ~isempty(header)
            settings.data = header2string(header);
            settings.Filename = fullfile(cfg.EEG_filepath,cfg.EEG_file);
        else
            settings.Filename = '';
        end
        settings.data =  header2string(header);
    end
    
    function settings = do_edit(settings)
        List = {'None','Preprocessing','Extra channels','Exclude channels','Reduce Length','Edit LOCS','Interpolate bad','Re-reference data'};
        [selection] = listdlg('PromptString','Select','Name','Edit','SelectionMode','single', ...
            'ListString',List,'CancelString','None','ListSize',[150 120]);
        if selection == 2
            settings = do_preprocessing(settings);
        elseif selection == 3
            settings = set_extra(settings);
        elseif selection == 4
            settings = exclude_channels(settings);
        elseif selection == 5
            settings = reduce_datalength(settings);
        elseif selection == 6
            settings = edit_locs(settings);
        elseif selection == 7
            settings = interpolate_bad(settings);
        elseif selection == 8
            settings = reference_data(settings);
        end
    end

    function settings = do_preprocessing(settings)
        if isempty(header)
            return
        end
        cfg.FILT = [];
        cfg.SCALE = [];
        cfg.RESAMPLE = [];
        cfg.BADELEC = [];
        cfg.MARK = [];
        cfg = lab_set_preprocessing(cfg,header,0,1);
        [data,header,cfg] = lab_preprocessing(data,header,cfg);
        settings.data =  header2string(header);
    end

    function settings = set_extra(settings)
        if isempty(header)
            return
        end
        cfg = lab_set_extrachannels(cfg,header,data);
        [data,header,cfg] = lab_extrachannels(data,header,cfg);
        settings.data =  header2string(header);
    end

    function settings = exclude_channels(settings)
        if isempty(header)
            return
        end
        cfg = lab_set_exclude(cfg,header);
        if isfield(cfg,'exclude') & ~isempty(cfg.exclude)
            include = setdiff(1:header.numchannels,cfg.exclude);
            [data,header] = lab_reduce_channels(data,header,include);
            settings.data =  header2string(header);
        end
    end

    function settings = reduce_datalength(settings)
        if isempty(header)
            return
        end
        [cfg,skipprocessing] = lab_set_reduce_datalength(cfg,header);
        if skipprocessing == 1
            return
        end
        [data,header,cfg] = lab_reduce_datalength(data,header,cfg);
        settings.data =  header2string(header);
    end

    function settings = edit_locs(settings)
        if isempty(header)
            return
        end
        if ~isfield(header,'locs') | isempty(header.locs)
            if ~isfield(header,'numdatachannels')
                header.numdatachannels = header.numchannels;
            end
            if ~isfield(header,'locs')
                header.locs = [];
            end
            header.locs = lab_load_locs(header.locs,cfg,header.numdatachannels);
        else
            settings2.LOCS = locs;
            [~,~,locs2] = lab_plot_locs(settings2,1,1,0,0);
            header.locs = locs2;
        end
        settings.data =  header2string(header);
    end

    function settings = interpolate_bad(settings)
        if isempty(header)
            return
        end
        [settings,skipprocessing] = lab_set_interpolate_bad(settings,header);
        if skipprocessing == 1
            return
        else
            header.badchans = settings.badchannels;
        end
        [data,header] = lab_interpolate_bad(data,header,settings.method,1);
        settings.data =  header2string(header);
    end

    function settings = reference_data(settings)
        if isempty(header)
            return
        end
        [cfg,skipprocessing] = lab_set_references(cfg,header);
        if skipprocessing == 1
            return
        end
        [data,header,cfg] = lab_reference_data(data,header,cfg);
        settings.data =  header2string(header);
    end

    
end

function Fout = header2string(header)
  if isempty(header)
      Fout = '';
      return
  end
  Fout = {};
  Fout{end+1,1} = ['Channels: ' num2str(header.numchannels)];
  if isfield(header,'numdatachannels')
      Fout{end+1,1} = ['Data channels: ' num2str(header.numdatachannels)];
  end
  if isfield(header,'numauxchannels')
      Fout{end+1,1} = ['Auxillary channels: ' num2str(header.numauxchannels)];
  end
  Fout{end+1,1} = ['Samplingrate: ' num2str(header.samplingrate)];
  Fout{end+1,1} = ['Length (timeframes): ' num2str(header.numtimeframes)];
  if isfield(header,'ref_chan') & isnumeric(header.ref_chan) & ~isempty(header.ref_chan)
      Fout{end+1,1} = ['Reference channel: ' num2str(header.ref_chan)];
  elseif isfield(header,'ref_chan') & ischar(header.ref_chan) & ~isempty(header.ref_chan)
      Fout{end+1,1} = ['Reference: ' header.ref_chan];
  else
      Fout{end+1,1} = ['Reference: unkown'];
  end
  if isfield(header,'ecg_ch') & ~isempty(header.ecg_ch)
      Fout{end+1,1} = ['ECG channel: ' num2str(header.ecg_ch)];
  end
  if isfield(header,'eog_ch') & ~isempty(header.eog_ch)
      Fout{end+1,1} = ['EOG channel: ' num2str(header.eog_ch)];
  end
  Fout = char(Fout);
end