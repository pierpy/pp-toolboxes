% Plot spectras 
%
% lab_evaluate_spectras
%
%     Input         Spectra.mat-file (lab_calculate_spectras)
%
% Written by F. Hatz 2013

function lab_evaluate_spectras

[settings,skipprocessing] = lab_set_evaluate_spectras;
if skipprocessing == 1
    return
else
    pause(0.2);
end

if isempty(settings.Spect)
    return
end

[~,calc.result_path] = lab_filename(settings.file);
calc.resultplots_path = fullfile(calc.result_path,'Plots');
calc.result_path = calc.resultplots_path;
warning off %#ok<WNOFF>
mkdir(calc.resultplots_path);
warning on %#ok<WNON>

if isfield(settings.Spect,'header')
    header = settings.Spect.header;
else
    header = [];
end
if ~isfield(header,'locs') & ~isempty(settings.Locs)
    header.locs = settings.Locs;
end
if isfield(header,'patient')
    patient = header.patient;
else
    patient = lab_subjectname(settings.file);
end
tmp = regexp(patient,'\d');
if length(tmp) == length(patient)
    patient = ['P_' patient];
end
clearvars tmp
if isfield(settings.Spect,'channels')
    channels = settings.Spect.channels;
end

Mappings = settings.mappings;
if isfield(settings.Spect,'SpectAllValid')
    Valid = [num2str(size(settings.Spect.SpectAll,1)) ' of ' num2str(size(settings.Spect.SpectAllValid,2))];
else
    Valid = '';
end
SpectAll = settings.Spect.SpectAll;
SpectAllMedian = settings.Spect.SpectAllMedian;
SpectAllMean = settings.Spect.SpectAllMean;
SpectAllF = settings.Spect.SpectAllF;

if strcmp(settings.source,'single')
    MaxSpect = size(SpectAll,3);
    cfg.CollectFFT.source = 'median';
else
    MaxSpect = 1;
end

% Define names for output-files
calc.output_fileBA = 'BackgroundActivity';
calc.output_fileQ = 'QualityCheck';
calc.output_filePF = 'MappingPeakFreq';
calc.output_fileMF = 'MappingMedianFreq';
calc.output_fileEPF = 'PeakFreq';
calc.output_fileEMF = 'MedianFreq';
calc.output_fileBP = 'MappingBandPowerMean';
calc.output_fileBPR = 'MappingBandPowerMeanRelative';
calc.output_fileBPM = 'MappingBandPowerMedian';
calc.output_fileBPMR = 'MappingBandPowerMedianRelative';
calc.output_fileEBP = 'BandPowerMean';
calc.output_fileEBPR = 'BandPowerMeanRelative';
calc.output_fileEBPM = 'BandPowerMedian';
calc.output_fileEBPMR = 'BandPowerMedianRelative';
calc.matfile = settings.file;
calc.mappingBA = settings.mappingBA;
calc.SpectAll = SpectAll;
calc.SpectAllF = SpectAllF;
cfg.CollectFFT = settings;
if isfield(cfg.CollectFFT,'spectralbandsI') & cfg.CollectFFT.spectralbandsI == true & isfield(header,'IFREQ') & isfield(header.IFREQ,'Bands')
    cfg.CollectFFT.freqs = cell2mat(header.IFREQ.Bands(:,2:5));
elseif ~isempty(cfg.CollectFFT.spectralbands)
    cfg.CollectFFT.freqs = cell2mat(cfg.CollectFFT.spectralbands(:,2:3));
    cfg.CollectFFT.freqs = [cfg.CollectFFT.freqs cfg.CollectFFT.freqs];
end
[cfg,result,calc,Mappings,skipprocessing] = lab_collect_spectraldata_header(cfg,calc,[],Mappings);
if skipprocessing == 1
    return
end

result.FBM = 'empty';
for Nspect = 1:MaxSpect
    if MaxSpect > 1
        calc.patient = [patient '_E' num2str(Nspect)];
        SpectAllMedian =permute(SpectAll(Nspect,:,:),[3 2 1]);
        SpectAllMean = permute(SpectAll(Nspect,:,:),[3 2 1]);
    else
        calc.patient = patient;
    end
    % Calculate peak frequency & median frequency
    [cfg,result] = lab_collect_spectraldata_pfmf(cfg,result,calc,SpectAllF,SpectAllMean,SpectAllMedian,cfg.CollectFFT.mappings,Valid,header);
    if isnumeric(SpectAllF)
        % Calculate BandPower Mapping
        [cfg,result] = lab_collect_spectraldata_bp(cfg,result,calc,SpectAll,SpectAllMean,SpectAllMedian,SpectAllF,Mappings);
    end
    % Collect relative power
    [calc,result] = lab_collect_spectraldata_relativepower(cfg,calc,result,SpectAllMean,SpectAllMedian);
end

disp ('    Plot results')
result = lab_collect_spectraldata_plot(cfg,calc,result,false);
if settings.writeresults == true
    disp ('    Write results')
    lab_collect_spectraldata_write(cfg,calc,result,cfg.CollectFFT.mappings);
end
