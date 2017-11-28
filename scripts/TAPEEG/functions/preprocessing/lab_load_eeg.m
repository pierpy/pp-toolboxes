% Load eeg/meg-file for manual processing
%
% Written by F. Hatz 2013 Neurology Basel

function lab_load_eeg

if ~exist('cfg','var')
    cfg = [];
end
data = [];
header = [];
activations = [];
headerICA = [];
Spectra = [];
Epochs = [];
Connect = [];

settings.Filename = '';
settings.data = '';
settings.locs = [];
settings.badchans = [];
settings.interpolated = [];
settings.events = false;
settings.headerICA = headerICA;
settings.Spectra = [];
settings.Epochs = [];
settings.Connect = [];

Prompt = cell(0,2);
Formats = [];
Options.CancelButton = 'off';
Options.ApplyButton = 'off';
Options.OkButton = 'off';
Options.Wait = 'off';
Options.WindowStyle = 'normal';
Options.Menubar = 'on';

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

Prompt(end+1,:) = {'Locs','locs'};
Formats(end+1,1).type = 'check';
Formats(end,1).callback = {@edit_locs,'locs','locs'};

Prompt(end+1,:) = {'Marker','events'};
Formats(end+1,1).type = 'check';
Formats(end,1).format = 'integer';
Formats(end,1).callback = {@plot_eeg};

Prompt(end+1,:) = {'Bad channels','badchans'};
Formats(end+1,1).type = 'check';
Formats(end,1).callback = {@edit_badchans,'badchans','badchans'};

Prompt(end+1,:) = {'Interpolated channels','interpolated'};
Formats(end+1,1).type = 'check';
Formats(end,1).callback = {@show_interpolated,[],'interpolated'};

Prompt(end+1,:) = {'ICA activations','headerICA'};
Formats(end+1,1).type = 'check';
Formats(end,1).callback = {@do_ICA,'@ALL','@ALL'};

Prompt(end+1,:) = {' ',''};
Formats(end+1,1).type = 'text';
Formats(end,1).span = [1 2];

Prompt(end+1,:) = {'Results',''};
Formats(end+1,1).type = 'text';
Formats(end,1).span = [1 2];

Prompt(end+1,:) = {'Spectral analysis','Spectra'};
Formats(end+1,1).type = 'check';
Formats(end,1).callback = {@spectral_analysis};

Prompt(end+1,:) = {'Epochs','Epochs'};
Formats(end+1,1).type = 'check';
Formats(end,1).callback = {@calculate_epochs};

Prompt(end+1,:) = {'Connectivity','Connect'};
Formats(end+1,1).type = 'check';
Formats(end,1).callback = {@calculate_connectivity};

[~,handles] = inputsdlg(Prompt,'EEG/MEG',Formats,settings,Options);
H.data = handles.ctrl(3);
H.locs = handles.ctrl(4);
H.events = handles.ctrl(5);
H.badchans = handles.ctrl(6);
H.interpolated = handles.ctrl(7);
H.spectra = handles.ctrl(11);
H.epochs = handles.ctrl(12);
H.connect = handles.ctrl(13);

m1 = uimenu(gcf,'Label','Channels');
uimenu(m1,'Label','Define extra channels','Callback',@set_extra);
uimenu(m1,'Label','Exclude channels','Callback',@exclude_channels);
uimenu(m1,'Label','Interpolate bad channels','Callback',@interpolate_bad);

m2 = uimenu(gcf,'Label','Locs');
uimenu(m2,'Label','Plot locations','Callback',@plot_locs);
uimenu(m2,'Label','Load locations','Callback',@load_locs);

m3 = uimenu(gcf,'Label','Preprocessing');
uimenu(m3,'Label','Preprocessing','Callback',@do_preprocessing);
uimenu(m3,'Label','Reduce length','Callback',@reduce_datalength);

m4 = uimenu(gcf,'Label','Postprocessing');
uimenu(m4,'Label','Re-reference','Callback',@reference_data);
uimenu(m4,'Label','Spectral analysis','Callback',@spectral_analysis);
uimenu(m4,'Label','Split in Epochs','Callback',@calculate_epochs);
uimenu(m4,'Label','Connectivity','Callback',@calculate_connectivity);

m5 = uimenu(gcf,'Label','Plot');
uimenu(m5,'Label','EEG','Callback',@plot_eeg);
uimenu(m5,'Label','ICA','Callback',@plot_ICA);
uimenu(m5,'Label','Connectivity','Callback',@plot_matrix);

% nested functions
    function settings2 = read_eeg(settings2)
        if exist(settings2.Filename,'file')
            Filename = settings2.Filename;
        else
            Filename = '';
        end
        [data,header,cfg] = lab_read_data(Filename,cfg);
        if ~isempty(data) & ~isempty(header)
            settings2.data = header2string(header);
            settings2.Filename = fullfile(cfg.EEG_filepath,cfg.EEG_file);
            if isfield(header,'locs')
                settings2.locs = header.locs;
            end
            if isfield(header,'badchans') & ~isempty(header.badchans)
                settings2.badchans = header.badchans;
            else
                settings2.badchans = [];
            end
            if isfield(header,'interpolated') & ~isempty(header.interpolated)
                settings2.interpolated = header.interpolated;
            else
                settings2.interpolated = [];
            end
            if isfield(header,'events') & ~isempty(header.events)
                settings2.events = true;
            else
                settings2.events = false;
            end
        else
            settings2.Filename = '';
        end
        update_header;
    end

    function settings2 = do_ICA(settings2)
        if ~isempty(activations)
            doICA = questdlg('Backtransform or Delete','ICA','Delete','Backtransform','Backtransform');
            if strcmp(doICA,'Delete')
                activations = [];
                headerICA = [];
                settings2.headerICA = [];
            elseif strcmp(doICA,'Backtransform')
                [cfg,skipprocessing] = lab_set_ICAback(cfg,header,1);
                if skipprocessing == 1
                    return
                end
                if isfield(cfg.ICABACK,'ACTIVATIONS') & ~isempty(cfg.ICABACK.ACTIVATIONS)
                    headerICA.badchans = cfg.ICABACK.ACTIVATIONS;
                elseif isfield(cfg.ICABACK,'BAD') & ~isempty(cfg.ICABACK.BAD)
                    if isfield(cfg.ICABACK.BAD,'ecg_ch') & cfg.ICABACK.BAD.ecg_ch > 0 & cfg.ICABACK.BAD.ecg_ch <= size(data,1)
                        cfg.ICABACK.BAD.ecg_chan = data(cfg.ICABACK.BAD.ecg_ch,:);
                    end
                    if isfield(cfg.ICABACK.BAD,'eog') & ~isempty(cfg.ICABACK.BAD.eog)
                        cfg.ICABACK.BAD = lab_calculate_eog(data,header,cfg.ICABACK.BAD);
                    end
                    if isfield(headerICA,'locs')
                        headerICA = rmfield(headerICA,'locs');
                    end
                    [headerICA.badchans,~,ICABACK.BAD] = lab_detect_bad(activations,headerICA,cfg.ICABACK.BAD,cfg,1);
                end
                if isfield(cfg.ICABACK,'dobacktransform') & cfg.ICABACK.dobacktransform == true
                    W = headerICA.W;
                    if max(headerICA.badchans) > 0
                        W(:,headerICA.badchans) = 0;
                    end
                    data(header.goodchans,:)=W*activations;
                end
                if isfield(cfg.ICABACK,'IS') & ~isempty(cfg.ICABACK.IS)
                    include = setdiff(1:size(activations,1),headerICA.badchans);
                    W = W(:,include);
                    headerICA.W = W;
                    [activations2,headerICA2] = lab_reduce_channels(activations,headerICA,include,true);
                    lab_ICA_dipolfit(activations2,headerICA2,cfg);
                end
            end
        else
            if isempty(data) | isempty(header)
                return
            end
            if ~isfield(header,'goodchans')
                header.goodchans = 1:header.numdatachannels;
            end
            if isnumeric(header.ref_chan)
                header.goodchans = setdiff(header.goodchans,header.ref_chan);
                header.goodchans = header.goodchans(:)';
            end
            [activations,headerICA,cfg] = lab_ICAstart(data,header,cfg);
            if ~isempty(headerICA)
                settings2.headerICA = headerICA;
            end
        end
    end
    
    function plot_eeg(H1,H2)
        if isempty(header)
            return
        end
        [data,header,cfg] = lab_plot_eeg(data,header,cfg);
        update_header;
    end

    function plot_ICA(H1,H2)
        if ~isempty(activations)
            [activations,headerICA,cfg] = lab_plot_eeg(activations,headerICA,cfg);
        end
    end

    function do_preprocessing(H1,H2)
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
        update_header;
    end

    function set_extra(H1,H2)
        if isempty(header)
            return
        end
        cfg = lab_set_extrachannels(cfg,header,data);
    end

    function exclude_channels(H1,H2)
        if isempty(header)
            return
        end
        cfg = lab_set_exclude(cfg,header);
        if isfield(cfg,'exclude') & ~isempty(cfg.exclude)
            include = setdiff(1:header.numchannels,cfg.exclude);
            [data,header] = lab_reduce_channels(data,header,include);
            update_header;
        end
    end

    function plot_locs(H1,H2)
        if isfield(header,'locs') & ~isempty(header.locs)
            settings.LOCS = header.locs;
            lab_plot_locs(settings,1,0,0);
        end
    end

    function reduce_datalength(H1,H2)
        if isempty(header)
            return
        end
        [cfg,skipprocessing] = lab_set_reduce_datalength(cfg,header);
        if skipprocessing == 1
            return
        end
        [data,header,cfg] = lab_reduce_datalength(data,header,cfg);
    end

    function load_locs(H1,H2)
        if ~isfield(header,'numdatachannels')
            header.numdatachannels = header.numchannels;
        end
        if ~isfield(header,'locs')
            header.locs = [];
        end
        header.locs = lab_load_locs(header.locs,cfg,header.numdatachannels);
        update_header;
    end

    function locs = edit_locs(locs)
        if isempty(locs)
            if ~isfield(header,'numdatachannels')
                header.numdatachannels = header.numchannels;
            end
            if ~isfield(header,'locs')
                header.locs = [];
            end
            locs = lab_load_locs(header.locs,cfg,header.numdatachannels);
            header.locs = locs;
        else
            settings2.LOCS = locs;
            [~,~,locs2] = lab_plot_locs(settings2,1,1,0,0);
            header.locs = locs2;
            locs = locs2;
        end
    end

    function badchans = edit_badchans(badchans)
        settings2.indexed = badchans;
        if isfield(header,'locs') & ~isempty(header.locs)
            settings2.LOCS = header.locs;
            badchans = lab_plot_locs(settings2,1,0,0);
            header.badchans = badchans;
        else
            Prompt2 = {'Bad Channels','indexed'};
            Formats2.type = 'edit';
            Formats2.format = 'vector';
            Formats2.limits = [-inf inf];
            Formats2.size = 300;
            settings2 = inputsdlg(Prompt2,'Bad channels',Formats2,settings2);
            badchans = settings2.indexed;
        end
    end

    function show_interpolated(interpolated)
        questdlg(['Interpolated channels: ' num2str(interpolated)],'Interpolated','OK','OK');
    end

    function interpolate_bad(H1,H2)
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
        update_header;
    end

    function reference_data(H1,H2)
        if isempty(header)
            return
        end
        [cfg,skipprocessing] = lab_set_references(cfg,header);
        if skipprocessing == 1
            return
        end
        [data,header,cfg] = lab_reference_data(data,header,cfg);
        update_header;
    end

    function spectral_analysis(H1,H2)
        if isempty(header)
            return
        end
        if ~isempty(Spectra)
            doSpectra = questdlg('Show result or Delete','ICA','Delete','Show result','Show result');
            if strcmp(doSpectra,'Delete')
                Spectra = [];
                set(H_spectra,'Value',0);
                uicontrol(H_spectra);
            else
                lab_show_spectraresults(Spectra,cfg.FFT.Mappings);
            end
        else
            [cfg,skipprocessing] = lab_set_calculate_spectras(cfg,header);
            pause(0.2);
            if skipprocessing == 1
                return
            end
            [Spectra,cfg] = lab_calculate_spectras(data,header,cfg);
            set(H.spectra,'Value',1);
            uicontrol(H.spectra);
        end
    end

    function calculate_epochs(H1,H2)
        if isempty(header)
            return
        end
        if isempty(Epochs)
            [cfg,skipprocessing] = lab_set_save_epochs(cfg,header);
            pause(0.2);
            if skipprocessing == 1
                return
            end
            [Epochs,cfg] = lab_save_epochs(data,header,cfg);
            update_header;
            set(H.epochs,'Value',1);
            uicontrol(H.epochs);
        else
            doEpochs = questdlg('Delete?','Epochs','No','No');
            if strcmp(doEpochs,'Delete')
                Epochs = [];
                set(H.epochs,'Value',0);
                uicontrol(H.epochs);
            end
        end
    end

    function calculate_connectivity(H1,H2)
        if isempty(header)
            return
        end
        if isempty(Connect)
            if ~isempty(Epochs)
                doEpochs = questdlg('Calculate Connectivity for input data or epochs?','Connectivity','Inputdata','Epochs','Epochs');
            else
                doEpochs = 'Inputdata';
            end
            if strcmp(doEpochs,'Epochs') & ~isempty(Epochs)
                headerE = header;
                headerE.numtimeframes = size(Epochs.data,2);
                [cfg,skipprocessing] = lab_set_calculate_connectivity(cfg,headerE,Epochs.data(:,:,1),1);
                pause(0.2)
                if skipprocessing == 1
                    return
                end
                for i = 1:size(Epochs.data,3)
                    headerE.badchans = Epochs.badchans{1,i};
                    headerE.goodchans = setdiff(1:headerE.numdatachannels,headerE.badchans);
                    headerE.badchans = headerE.badchans(:)';
                    headerE.goodchans = headerE.goodchans(:)';
                    headerE.timestamp = Epochs.events.POS(1,i);
                    [Connect(i),cfg] = lab_connect(Epochs.data(:,:,i),headerE,cfg);
                end
            else
                [cfg,skipprocessing] = lab_set_calculate_connectivity(cfg,header,data,1);
                pause(0.2)
                if skipprocessing == 1
                    return
                end
                [Connect,cfg] = lab_connect(data,header,cfg,1);
                update_header;
                set(H.connect,'Value',1);
                uicontrol(H.connect);
            end
        else
            doConnect = questdlg('Delete?','Matrix','Delete','No','No');
            if strcmp(doConnect,'Delete')
                Connect = [];
                set(H.connect,'Value',0);
                uicontrol(H.connect);
            end
        end
    end

    function plot_matrix(H1,H2)
        if ~isempty(Connect);
            ResultS = Connect(1);
            content = fieldnames(ResultS);
            tmp = 0;
            for i = 1:size(content,1)
                if strcmp(content{i}(1),'F')
                    tmp = tmp+1;
                end
            end
            if tmp == size(content,1)
                selectionF = listdlg('ListString',content,'SelectionMode','single');
                if isempty(selectionF)
                    return
                end
                ResultS = Connect(1).(content{selectionF});
                content = fieldnames(ResultS);
            end
            tmp = [];
            for i = 1:size(content,1)
                if isnumeric(ResultS.(content{i})) & size(ResultS.(content{i}),1) == size(ResultS.(content{i}),2)
                    tmp = [tmp i];
                end
            end
            if isempty(tmp)
                return
            end
            selection = listdlg('ListString',content,'SelectionMode','single');
            if isempty(selection)
                return
            end
            selection = tmp(selection);
            
            Result = [];
            for i = 1:length(Connect)
                if exist('selectionF','var')
                    Resulttmp = Connect(i).(content{selectionF});
                else
                    Resulttmp = Connect(i);
                end
                Result = cat(3,Result,Resulttmp.(content{selection}));
            end
            Result = mean(Result,3);

            if isfield(header,'locs') & ~isempty(header.locs)
                doConnect = questdlg('Matrix or Connectivity-Plot','Matrix','Connectivity-Plot','Matrix','Matrix');
                if strcmp(doConnect,'Connectivity-Plot')
                    doConnect = questdlg('Signal or Source space','Connectivity-Plotx','Source','Signal','Signal');
                    if strcmp(doConnect,'Signal')
                        lab_plot_elec(Result,header.locs)
                    else
                        lab_plot_IS(Result);
                    end
                    return
                end
            end
            lab_plot_matrix(Result);
        end 
    end

    function update_header
        Fout = header2string(header);
        set(H.data,'String',Fout);
        uicontrol(H.data);
        clearvars Fout
        if isfield(header,'locs') & ~isempty(header.locs)
            set(H.locs,'Value',1);
        else
            set(H.locs,'Value',0);
        end
        uicontrol(H.locs);
        if isfield(header,'events') & ~isempty(header.events)
            set(H.events,'Value',1);
        else
            set(H.events,'Value',0);
        end
        uicontrol(H.events);
        if isfield(header,'badchans') & ~isempty(header.badchans)
            set(H.badchans,'Value',1,'UserData',header.badchans);
        else
            set(H.badchans,'Value',0,'UserData',[]);
        end
        uicontrol(H.badchans);
        if isfield(header,'interpolated') & ~isempty(header.interpolated)
            set(H.interpolated,'Value',1,'UserData',header.interpolated);
        else
            set(H.interpolated,'Value',0,'UserData',[]);
        end
        uicontrol(H.interpolated);
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
  Fout{end+1,1} = ['Length (seconds): ' num2str(header.numtimeframes/header.samplingrate,2)];
  if isnumeric(header.ref_chan) & ~isempty(header.ref_chan)
      Fout{end+1,1} = ['Reference channel: ' num2str(header.ref_chan)];
  elseif ischar(header.ref_chan) & ~isempty(header.ref_chan)
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