% Calculate new references and store results
%
% [epochs,cfg] = lab_save_epochs(data,header,cfg)
%
% data    = matrix (chans x timeframes)
% header  = output of lab_read_data
% cfg     = structure with config (optional)
%
% written by F. Hatz 2012

function [epochs,cfg] = lab_save_epochs(data,header,cfg)

if ~exist('header','var') | ~isfield(header,'samplingrate')
    header = lab_create_header(data);
end
if ~exist('cfg','var') || ~isfield(cfg,'settings_path')
    cfg.settings_path = pwd;
end
if ~isfield(cfg,'Output_file')
    cfg.Output_file = header.EEG_file;
    cfg.Output_filepath = header.EEG_filepath;
end
[~,~,~,cfg.Output_fileS] = lab_filename(cfg.Output_file);

if isfield(cfg,'SKIP') & isfield(cfg.SKIP,'EPOCH') & cfg.SKIP.EPOCH == true;
    epochs = [];
    return
end

disp('Save epochs')
if ~isfield(cfg,'EPOCH') | isempty(cfg.EPOCH)
    [cfg,skipprocessing] = lab_set_save_epochs(cfg,header);
    if skipprocessing == 1
        epochs = [];
        return
    end
end

cfg.EPOCH = lab_markerselect(cfg.EPOCH,cfg,header);

% Interpolate bad channels for analysis
if ~isfield(header,'goodchans')
    header.goodchans = 1:header.numdatachannels;
    header.badchans = [];
end
if cfg.EPOCH.interpolatebad == 1
    [data,header] = lab_interpolate_bad(data,header);
end

% Calculate new reference
if isfield(cfg.EPOCH,'eegsource') & ~strcmp(cfg.EPOCH.eegsource,'input')
    if ~isfield(header,'ref_chan') | (~isnumeric(header.ref_chan) & ~strcmp(header.ref_chan,'none'))
        epochs = [];
        return
    end
    if strcmp(cfg.EPOCH.eegsource,'montage') & isfield(cfg.EPOCH,'montage') & ~isempty(cfg.EPOCH.montage)
        [data,header,cfg.EPOCH] = lab_references(data,header,cfg.EPOCH.montage,cfg.EPOCH);
    else
        [data,header,cfg.EPOCH] = lab_references(data,header,cfg.EPOCH.eegsource,cfg.EPOCH);
    end
end

% Find epochs
if size(data,2) > cfg.EPOCH.length
    settings = cfg.EPOCH;
    [epochs,settings] = lab_select_epochs(data,header,settings);
    if isempty(epochs)
        disp('   Creating epochs failed')
        return
    end
    disp(['    found ' num2str(size(epochs.data,3)) ' epochs with ' num2str(epochs.percentgood) ' percent good channels'])
    cfg.EPOCH = settings;
    clearvars settings
    if cfg.EPOCH.interpolate3D == 1
        disp('   3D-Interpolation of bad channels (epochs)')
        for i = 1:size(epochs.data,3)
            headerepochs = header;
            if isfield(cfg.EPOCH,'replacebad') & cfg.EPOCH.replacebad == true
                headerepochs.badchans = union(header.badchans,epochs.badchans{1,i});
                headerepochs.badchans = headerepochs.badchans(:)';
                headerepochs.goodchans = setdiff(headerepochs.goodchans,headerepochs.badchans);
                headerepochs.goodchans = headerepochs.goodchans(:)';
            end
            [epochs.data(:,:,i),headerepochs] = lab_interpolate_bad(epochs.data(:,:,i),headerepochs,'3D');
            if isfield(headerepochs,'interpolated')
                epochs.interpolated{1,i} = headerepochs.interpolated;
            end
            epochs.badchans{1,i} = [];
        end
    else
        for i = 1:size(epochs.data,3)
            if isfield(header,'badchans') & (header.badchans)
                epochs.badchans{1,i} = header.badchans;
            end
            if isfield(header,'interpolated') & (header.interpolated)
                epochs.interpolated{1,i} = header.interpolated;
            end
        end
    end
    disp('   Write epochs to files')
    [epochs,cfg] = lab_write_epochs(epochs,header,cfg);
    header = lab_mix_markers(header,epochs);
    
    % -- Do calculations on epochs --
    headerE = header;
    headerE.numtimeframes = size(epochs.data,2);
    if isfield(headerE,'bad')
        headerE = rmfield(headerE,'bad');
    end
    if ~isfield(headerE,'cov')
        headerE.cov = cov(data(1:headerE.numdatachannels,:)');
    end
    Output_file = cfg.Output_file;
    Output_fileS = cfg.Output_fileS;
    Output_filepath = cfg.Output_filepath;
    cfg.Output_filepath = fullfile(cfg.Output_filepath,cfg.EPOCH.folder);
    
    if isfield(cfg,'lastsegment')
        lastsegment = cfg.lastsegment;
    else
        lastsegment = true;
    end
    if isfield(cfg,'lastfilefolder')
        lastfilefolder = cfg.lastfilefolder;
    else
        lastfilefolder = true;
    end
    if isfield(cfg,'firstfilefolder')
        firstfilefolder = cfg.firstfilefolder;
    else
        firstfilefolder = true;
    end
    cfg.SKIP.EPOCH = true;
    SKIP = cfg.SKIP;
    for i = 1:size(epochs.data,3)
        if isfield(epochs,'goodsum') & length(epochs.goodsum) >= i
            headerE.quality = epochs.goodsum(i);
        end
        if lastsegment == true & i == size(epochs.data,3)
            cfg.lastsegment = true;
        else
            cfg.lastsegment = false;
        end
        if lastfilefolder == true & i == size(epochs.data,3)
            cfg.lastfilefolder = true;
        else
            cfg.lastfilefolder = false;
        end
        if firstfilefolder == true & i == 1
            cfg.firstfilefolder = true;
        else
            cfg.firstfilefolder = false;
        end
        cfg.Output_file = [Output_fileS '_E' num2str(i) '.sef'];
        cfg.Output_fileS = [Output_fileS '_E' num2str(i)];
        if isfield(cfg.EPOCH,'replacebad') & cfg.EPOCH.replacebad == true
            headerE.badchans = epochs.badchans{1,i};
            headerE.badchans = headerE.badchans(:)';
            headerE.goodchans = setdiff(1:headerE.numdatachannels,headerE.badchans);
            headerE.goodchans = headerE.goodchans(:)';
        end
        cfg.SKIP = SKIP;
        cfg = lab_postprocessing(epochs.data(:,:,i),headerE,cfg,true);
    end
    cfg.lastsegment = lastsegment;
    cfg.lastfilefolder = lastfilefolder;
    cfg.firstfilefolder = firstfilefolder;
    cfg.Output_file = Output_file;
    cfg.Output_fileS = Output_fileS;
    cfg.Output_filepath = Output_filepath;
end

%--------------------------------------------------------------------------
% Write verbose file (*.vrb)
%--------------------------------------------------------------------------
cd(fullfile(cfg.Output_filepath,cfg.EPOCH.folder));
fid=fopen([cfg.Output_fileS '.vrb'],'w');
fprintf(fid,'Save Epochs\n');
fprintf(fid,['EEG-file: ' cfg.Output_file]);
fprintf(fid,'\n\n');
fprintf(fid,['Length of exported epochs (timeframes): ' num2str(cfg.EPOCH.length)]);
fprintf(fid,'\n');
if isfield(cfg.EPOCH,'BAD')
    fprintf(fid,'Settings for detecting bad channels:\n');
    fprintf(fid,['  Bad Spectra 50Hz (Limit ' num2str(round(cfg.EPOCH.BAD.freqlim50)) 'percent)\n']);
    fprintf(fid,['  Bad Spectra 60Hz (Limit ' num2str(round(cfg.EPOCH.BAD.freqlim60)) 'percent)\n']);
    fprintf(fid,['  Bad Spectra '  num2str(cfg.EPOCH.BAD.spectshigh(1,1)) '-' num2str(cfg.EPOCH.BAD.spectshigh(1,2)) ...
        'Hz (Limit ' num2str(round(cfg.EPOCH.BAD.freqlimhigh)) 'percent)\n']);
    fprintf(fid,['  Bad Spectra '  num2str(cfg.EPOCH.BAD.spectslow(1,1)) '-' num2str(cfg.EPOCH.BAD.spectslow(1,2)) ...
        'Hz (Limit ' num2str(round(cfg.EPOCH.BAD.freqlimlow)) 'percent)\n']);
    fprintf(fid,['  Bad variance (zValue: '  num2str(cfg.EPOCH.BAD.zvaluevars) ')\n']);
    fprintf(fid,['  Bad hurst (zValue: '  num2str(cfg.EPOCH.BAD.zvaluehurst) ')\n']);
    fprintf(fid,['  Bad topo correlation (zValue: '  num2str(cfg.EPOCH.BAD.zvaluecorr) ')\n']);
    if isfield(cfg.EPOCH.BAD,'PEAK2MIN') & ~isempty(cfg.EPOCH.BAD.PEAK2MIN)
        if strcmp(cfg.EPOCH.BAD.PEAK2MIN.mode,'threshold')
            fprintf(fid,['Peak2min (' num2str(cfg.EPOCH.BAD.PEAK2MIN.lowfreqpeak) '-' ...
                num2str(cfg.EPOCH.BAD.PEAK2MIN.highfreqpeak) 'Hz - ' ...
                cfg.EPOCH.BAD.PEAK2MIN.mode ...
                ' - factor:' num2str(cfg.EPOCH.BAD.PEAK2MIN.factor) ' - min:' ...
                num2str(cfg.EPOCH.BAD.PEAK2MIN.threshold) '):\n']);
        else
            fprintf(fid,['Peak2min (' num2str(cfg.EPOCH.BAD.PEAK2MIN.lowfreqpeak) '-' ...
                num2str(cfg.EPOCH.BAD.PEAK2MIN.highfreqpeak) 'Hz - ' ...
                cfg.EPOCH.BAD.PEAK2MIN.mode ...
                ' - factor:' num2str(cfg.EPOCH.BAD.PEAK2MIN.factor) '):\n']);
        end
        fprintf(fid,'\n');
    end
end
fprintf(fid,['Percent of good channels per epoch: ' num2str(epochs.percentgood) '%']);
fprintf(fid,'\n');
if cfg.EPOCH.minimalpart <= 1
    fprintf(fid,['Fraction of summed epochs on eeg length >= ' num2str(cfg.EPOCH.minimalpart*100) '%']);
else
    fprintf(fid,['Number of epochs: ' num2str(cfg.EPOCH.minimalpart)]);
end
fprintf(fid,'\n');
if cfg.EPOCH.interpolatebad == 1
    fprintf(fid,'Interpolate bad channels: enabled');
else
    fprintf(fid,'Interpolate bad channels: disabled');
end
fprintf(fid,'\n');
if ~isempty(cfg.EPOCH.markerexclude)
    fprintf(fid,['Markers excluded: ' sprintf('%s|',cfg.EPOCH.markerexclude{:})]);
    fprintf(fid,'\n');
end
if ~isempty(cfg.EPOCH.markerinclude)
    fprintf(fid,['Markers included: ' sprintf('%s|',cfg.EPOCH.markerinclude{:})]);
    fprintf(fid,'\n');
end
if isfield(cfg.EPOCH,'scaletxtnum')
    fprintf(fid,['Scalefactor used for txt-export: ' num2str(cfg.EPOCH.scaletxtnum)]);
elseif strcmp(cfg.EPOCH.scaletxt,'none')
    fprintf(fid,'Scalefactor used for txt-export: none');
end
fclose(fid);

if ~isempty(cfg.EPOCH.markerexclude) & strcmp(cfg.EPOCH.markerexclude{1,1},'all')
    cfg.EPOCH.markerexclude = cellstr('all');
end
if ~isempty(cfg.EPOCH.markerinclude) & strcmp(cfg.EPOCH.markerinclude{1,1},'all')
    cfg.EPOCH.markerinclude = cellstr('all');
end

return