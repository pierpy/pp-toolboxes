% Script to run postprocessings, called by lab_tapeeg
%
% Written by F. Hatz 2012/2014

function cfg = lab_postprocessing(data,header,cfg,noskipreset)

if ~exist('noskipreset','var') | isempty(noskipreset)
    noskipreset = false;
end
if noskipreset == false
    disp ('-Postprocessing-')
    cfg.SKIP = [];
elseif ~exist('cfg','var')
    cfg = [];
end
if ~isfield(cfg,'Output_file')
    if isfield(cfg,'EEG_file')
        cfg.Output_file  = cfg.EEG_file;
        cfg.Output_filepath  = cfg.EEG_filepath;
    elseif isfield(header,'EEG_file')
        cfg.Output_file  = header.EEG_file;
        cfg.Output_filepath  = header.EEG_filepath;
    end
end
if isfield(cfg,'Output_file')
    [~,~,~,cfg.Output_fileS] = lab_filename(cfg.Output_file);
end
    
% Write result file
if isfield(cfg,'RES') & ~isempty(cfg.RES)
    cfg = lab_write_result(data,header,cfg);
end

% Write edf file
if isfield(cfg,'EDF') & ~isempty(cfg.EDF)
    cfg.EDF_file = [cfg.Output_fileS '.edf'];
    cfg.EDF_filepath = cfg.Output_filepath;
    [~,~,cfg] = lab_export2edf(data,header,cfg);
end

% Reduce channels
if isfield(cfg,'exclude') & ~isempty(cfg.exclude) & (~isfield(cfg.SKIP,'EXCLUDE') | cfg.SKIP.EXCLUDE == false)
    disp (['  Exclude predefined ' num2str(length(cfg.exclude)) ' channels'])
    [data,header] = lab_reduce_channels(data,header,setdiff(1:header.numchannels,cfg.exclude));
    cfg.SKIP.EXCLUDE = true;
end

% Re-reference data
if isfield(cfg,'REF') & ~isempty(cfg.REF)
    [~,~,cfg] = lab_reference_data(data,header,cfg);
end

% Epoch data
if isfield(cfg,'EPOCH') & ~isempty(cfg.EPOCH)
    [~,cfg] = lab_save_epochs(data,header,cfg);
end

% Average
if isfield(cfg,'AVG') & ~isempty(cfg.AVG)
    [~,cfg] = lab_average(data,header,cfg);
end

% Inverse solution
if isfield(cfg,'IS') & isfield(cfg.IS,'eegsource')
    if noskipreset == false
        if exist('lab_inversesolution') %#ok<EXIST>
            [~,cfg] = lab_inversesolution(data,header,cfg);
        else
            disp('Skip - no implemented yet')
        end
    elseif strcmp(cfg.IS.eegsource,header.ref_chan) | strcmp(cfg.IS.eegsource,'input')
        eegsource = cfg.IS.eegsource;
        cfg.IS.eegsource = 'input';
        if exist('lab_inversesolution') %#ok<EXIST>
            [~,cfg] = lab_inversesolution(data,header,cfg);
        else
            disp('Skip - no implemented yet')
        end
        cfg.IS.eegsource = eegsource;
        clearvars eegsource
    end
end

% Inverse solution fft
if isfield(cfg,'ISFFT') & isfield(cfg.ISFFT,'eegsource') & ~strcmp(header.ref_chan,'source')
    if noskipreset == false 
        if exist('lab_inversesolution_fft') %#ok<EXIST>
            [~,cfg] = lab_inversesolution_fft(data,header,cfg);
        else
            disp('Skip - no implemented yet')
        end
        cfg.SKIP.ISFFT = true;
    elseif strcmp(cfg.ISFFT.eegsource,header.ref_chan) | strcmp(cfg.ISFFT.eegsource,'input')
        eegsource = cfg.ISFFT.eegsource;
        cfg.ISFFT.eegsource = 'input';
        if exist('lab_inversesolution_fft') %#ok<EXIST>
            [~,cfg] = lab_inversesolution_fft(data,header,cfg);
        else
            disp('Skip - no implemented yet')
        end
        cfg.ISFFT.eegsource = eegsource;
        clearvars eegsource
        cfg.SKIP.ISFFT = true;
    end
end

% Spectras
if isfield(cfg,'FFT') & isfield(cfg.FFT,'eegsource')
    if noskipreset == false
        [~,cfg] = lab_calculate_spectras(data,header,cfg);
        cfg.SKIP.FFT = true;
    elseif strcmp(cfg.FFT.eegsource,header.ref_chan) | strcmp(cfg.FFT.eegsource,'input')
        eegsource = cfg.FFT.eegsource;
        cfg.FFT.eegsource = 'input';
        [~,cfg] = lab_calculate_spectras(data,header,cfg);
        cfg.FFT.eegsource = eegsource;
        clearvars eegsource
        cfg.SKIP.FFT = true;
    end
end

% Microstates
if isfield(cfg,'MICROST') & isfield(cfg.MICROST,'eegsource')
    if noskipreset == false
        if exist('lab_microstates') %#ok<EXIST>
            [~,~,cfg] = lab_microstates(data,header,cfg);
        else
            disp('Skip - no implemented yet')
        end
        cfg.SKIP.MICROST = true;
    elseif strcmp(cfg.MICROST.eegsource,header.ref_chan) | strcmp(cfg.MICROST.eegsource,'input')
        eegsource = cfg.MICROST.eegsource;
        cfg.MICROST.eegsource = {'input'};
        if exist('lab_microstates') %#ok<EXIST>
            [data,header,cfg] = lab_microstates(data,header,cfg);
        else
            disp('Skip - no implemented yet')
        end
        cfg.MICROST.eegsource = eegsource;
        clearvars eegsource
        cfg.SKIP.MICROST = true;
    end
end

% Connectivity
if isfield(cfg,'CONNECT') & isfield(cfg.CONNECT,'eegsource')
    if noskipreset == false 
        if exist('lab_calculate_connectivity') %#ok<EXIST>
            [~,cfg] = lab_calculate_connectivity(data,header,cfg);
        else
            disp('Skip - no implemented yet')
        end
        cfg.SKIP.CONNECT = true;
    elseif strcmp(cfg.CONNECT.eegsource,header.ref_chan) | strcmp(cfg.CONNECT.eegsource,'input')
        eegsource = cfg.CONNECT.eegsource;
        cfg.CONNECT.eegsource = 'input';
        if exist('lab_calculate_connectivity') %#ok<EXIST>
            [~,cfg] = lab_calculate_connectivity(data,header,cfg);
        else
            disp('Skip - no implemented yet')
        end
        cfg.CONNECT.eegsource = eegsource;
        clearvars eegsource
        cfg.SKIP.CONNECT = true;
    end
end

% Distance matrices
if isfield(cfg,'DISTANCE') & cfg.DISTANCE == true
    [~,cfg] = lab_calculate_distance(data,header,cfg);
    cfg.SKIP.DISTANCE = true;
end

% Forward solution
if isfield(cfg,'FWS') & ~isempty(cfg.FWS) & noskipreset == false
    if exist('lab_forwardsolution') %#ok<EXIST>
        [~,cfg] = lab_forwardsolution(data,header,cfg);
    else
        disp('Skip - no implemented yet')
    end
    cfg.SKIP.FWS = true;
end

% Coregistration
if isfield(cfg,'COREG') & ~isempty(cfg.COREG) & noskipreset == false
    if exist('lab_coregistration') %#ok<EXIST>
        [~,cfg] = lab_coregistration(data,header,cfg);
    else
        disp('Skip - no implemented yet')
    end
    cfg.SKIP.COREG = true;
end

end
