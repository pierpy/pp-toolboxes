% Collect results of frequency analysis (signal and inverse solution)
%
% cfg = lab_collect_spectraldata(cfg)
%
% cfg = structure with config (optional)
%
% written by F. Hatz 2013

function cfg = lab_collect_spectraldata(cfg)

disp('Collect spectral data')
if ~exist('cfg','var')
    cfg = [];
    skipselection = false;
else
    skipselection = true;
end
Files = {};

% Ask for settings
if ~isfield(cfg,'CollectFFT') | isempty(cfg.CollectFFT)
    [cfg,Files,skipprocessing] = lab_set_collect_spectraldata(cfg);
    if skipprocessing == 1
        return
    end
elseif skipselection == true
    % disable interactive windows for collecting in automatic mode
    cfg.CollectFFT.correctpf = false;
end

if isempty(Files)
    Files = lab_collect_spectraldata_search(cfg,skipselection);
end
if isempty(Files)
    return
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
calc.output_fileEBE = 'BandEntropyMean';
calc.output_fileETE = 'TotalEntropyMean';
calc.output_fileEBEM = 'BandEntropyMedian';
calc.output_fileETEM = 'TotalEntropyMedian';
calc.output_fileBE = 'MappingBandEntropyMean';
calc.output_fileTE = 'MappingTotalEntropyMean';
calc.output_fileBEM = 'MappingBandEntropyMedian';
calc.output_fileTEM = 'MappingTotalEntropyMedian';

% Do calculation for FFT-data
if ~isfield(cfg.CollectFFT,'folder')
    cfg.CollectFFT.folder = 'FrequencyAnalysis';
end
calc.result_path = fullfile(cfg.CollectFFT.searchfolder,cfg.CollectFFT.folder);
calc.resultplots_path = fullfile(calc.result_path,'Plots');
warning off %#ok<WNOFF>
mkdir(calc.result_path);
mkdir(calc.resultplots_path);
warning on %#ok<WNON>

% Set results first empty
result.BA = 'empty';
result.PF = 'empty';
result.MF = 'empty';
result.EPF = 'empty';
result.EMF = 'empty';
result.BP = 'empty';
result.BPR = 'empty';
result.BPM = 'empty';
result.BPMR = 'empty';
result.EBP = 'empty';
result.EBPR = 'empty';
result.EBPM = 'empty';
result.EBPMR = 'empty';
result.FBM = 'empty';
result.BPlog = 'empty';
result.BPRlog = 'empty';
result.BPMlog = 'empty';
result.BPMRlog = 'empty';
result.EBPlog = 'empty';
result.EBPRlog = 'empty';
result.EBPMlog = 'empty';
result.EBPMRlog = 'empty';
result.FBMpat = 'empty';
result.BE = 'empty';
result.TE = 'empty';
result.BEM = 'empty';
result.TEM = 'empty';
result.EBE = 'empty';
result.ETE = 'empty';
result.EBEM = 'empty';
result.ETEM = 'empty';
result.Q = 'empty';

% Prepare
calc.matfile = Files{1};

% Set electrodes for background activity
calc.mappingBA = cfg.CollectFFT.mappingBA;

% Calculate headers
[cfg,result,calc,cfg.CollectFFT.mappings,skipprocessing] = lab_collect_spectraldata_header(cfg,calc,result,cfg.CollectFFT.mappings);
if skipprocessing == 1
    return
end

for filenr = 1:length(Files);
    calc.matfile = Files{filenr};
    matfileLog = lab_filename(calc.matfile);
    disp (['    Collect ' matfileLog])
    
    % Load Spectras
    MAT = load(calc.matfile);
    if (~isfield(MAT,'header') | ~isfield(MAT.header,'locs')) & ~isempty(cfg.CollectFFT.Locs)
        MAT.header.locs = cfg.CollectFFT.Locs;
    end
    if isfield(MAT,'SpectAll') & isfield(MAT,'SpectAllF') & size(MAT.SpectAllF,2) == cfg.numfreqbins
        if ~isreal(MAT.SpectAll)
            MAT.SpectAll = abs(MAT.SpectAll);
        end
        if isfield(MAT,'SpectAllMean') & ~isempty(MAT.SpectAllMean) & ~isreal(MAT.SpectAllMean)
            MAT.SpectAllMean = abs(MAT.SpectAllMean);
        end
        if isfield(MAT,'SpectAllMedian') & ~isempty(MAT.SpectAllMedian) & ~isreal(MAT.SpectAllMedian)
            MAT.SpectAllMedian = abs(MAT.SpectAllMedian);
        end
        if isfield(MAT,'SpectAllCenter') & ~isempty(MAT.SpectAllCenter) & ~isreal(MAT.SpectAllCenter)
            MAT.SpectAllCenter = abs(MAT.SpectAllCenter);
        end
        
        if isfield(MAT,'header') & isfield(MAT.header,'patient')
            calc.patient = MAT.header.patient;
        end
        if ~isempty(cfg.CollectFFT.subjectname)
            [calc.patient,cfg.CollectFFT] = lab_subjectname(calc.matfile,cfg.CollectFFT);
        end
        tmp = regexp(calc.patient,'\d');
        if length(tmp) == length(calc.patient)
            calc.patient = ['P_' calc.patient];
        end
        clearvars tmp
        
        % Get individual frequency bands
        if isfield(cfg.CollectFFT,'spectralbandsI') & cfg.CollectFFT.spectralbandsI == true & isfield(MAT.header,'IFREQ') & isfield(MAT.header.IFREQ,'Bands')
            cfg.CollectFFT.freqs = cell2mat(MAT.header.IFREQ.Bands(:,2:5));
            if isnumeric(cfg.CollectFFT.spectralbands) & ~isempty(cfg.CollectFFT.spectralbands)
                cfg.CollectFFT.freqs = cfg.CollectFFT.freqs(cfg.CollectFFT.spectralbands,:);
            end
        elseif ~isempty(cfg.CollectFFT.spectralbands)
            cfg.CollectFFT.freqs = cell2mat(cfg.CollectFFT.spectralbands(:,2:3));
            cfg.CollectFFT.freqs = [cfg.CollectFFT.freqs cfg.CollectFFT.freqs];
        end
        tmp = cfg.CollectFFT.freqs(:,1:2);
        tmp(tmp>max(MAT.SpectAllF)) = max(MAT.SpectAllF);
        tmp(tmp<min(MAT.SpectAllF)) = min(MAT.SpectAllF);
        cfg.CollectFFT.freqs(:,1:2) = tmp;
        tmp = tmp(:,2) - tmp(:,1);
        tmp = find(tmp>0);
        if isempty(tmp)
            cfg.CollectFFT.freqs = [];
        else
            cfg.CollectFFT.freqs = cfg.CollectFFT.freqs(tmp,:);
        end
        clearvars tmp
        
        % Set number of valid epochs
        Valid = [num2str(size(MAT.SpectAll,1)) ' of ' num2str(size(MAT.SpectAllValid,2))];
        
        % Calculate peak frequency & median frequency
        [cfg,result] = lab_collect_spectraldata_pfmf(cfg,result,calc,MAT.SpectAllF,MAT.SpectAllMean,MAT.SpectAllMedian,cfg.CollectFFT.mappings,Valid,MAT.header);
        
        % Calculate BandPower Mapping
        if ~isempty(cfg.CollectFFT.freqs)
            [cfg,result] = lab_collect_spectraldata_bp(cfg,result,calc,MAT.SpectAll,MAT.SpectAllMean,MAT.SpectAllMedian,MAT.SpectAllF,cfg.CollectFFT.mappings);
        end
        
        % Collect relative power
        [calc,result] = lab_collect_spectraldata_relativepower(cfg,calc,result,MAT.SpectAllMean,MAT.SpectAllMedian);
        
        % Collect spectral entropy
        [cfg,result] = lab_collect_spectralentropy_bp(cfg,result,calc,MAT.SpectAll,MAT.SpectAllMean,MAT.SpectAllMedian,MAT.SpectAllF,cfg.CollectFFT.mappings);
        
        clearvars MAT
    else
        disp(['Skip ' calc.matfile ' -wrong input data'])
    end
end
% Write results
disp ('    Plot results')
result = lab_collect_spectraldata_plot(cfg,calc,result);
disp ('    Write results')
cfg = lab_collect_spectraldata_write(cfg,calc,result,cfg.CollectFFT.mappings);

cd(cfg.CollectFFT.searchfolder);

return