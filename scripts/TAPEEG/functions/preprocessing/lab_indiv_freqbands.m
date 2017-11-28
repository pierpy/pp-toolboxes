% Helper script for finding individual frequency bands
%
% [data,header,cfg] = lab_indiv_freqbands(data,header,cfg)
%
% Define individual frequency bands according to Moretti et al 2004
% - Peak frequency (PF) is found by first detection median frequency in given
%   frequency range, follow by detection of the peak at Median frequency +/- 1 Hz
% - Theta/alpha transition frequency (TF) is found by detection of the
%   lowest amplitude left to PF, higher as lowest frequency in detection range - 1 Hz
%
% written by F. Hatz 2016

function [data,header,cfg] = lab_indiv_freqbands(data,header,cfg)

if ~exist('cfg','var')
    cfg = [];
end
if ~exist('header','var') | isempty(header)
    disp('   Abort: header-information missing')
    return
end
if ~isfield(cfg,'IFREQ') | isempty(cfg.IFREQ)
    [cfg,skipprocessing] = lab_set_indiv_freqbands(cfg,header);
    if skipprocessing == 1
        return
    end
end

patient = '';
if isfield(header,'patient')
    patient = regexprep(header.patient,' ','_');
elseif isfield(cfg,'patient')
    patient = regexprep(cfg.patient,' ','_');
end

if isfield(header,'EEG_filepath') & ~isempty(patient)
    IFREQ_file = fullfile(header.EEG_filepath,[patient '.ifreq']);
    if exist(IFREQ_file,'file')
        try %#ok<TRYNC>
            load(IFREQ_file,'-mat');
            if exist('IFREQ','var')
                disp('   Individual frequency bands: read frequency bands from patient backup')
                header.IFREQ = IFREQ; %#ok<NODEF>
                return
            end
        end
    end
elseif isfield(header,'EEG_filepath')
    [~,~,~,EEG_fileS] = lab_filename(header.EEG_file);
    IFREQ_file = fullfile(header.EEG_filepath,[EEG_fileS '.ifreq']);
else
    IFREQ_file = '';
end

disp('   Calculate individual frequency bands')

if ~isfield(header,'numdatachannels')
    header.numdatachannels = size(data,1);
end

if ~isfield(cfg.IFREQ,'channels') | isempty(cfg.IFREQ.channels)
    Nact = 1:header.numdatachannels;
else
    Nact = cfg.IFREQ.channels;
end
winsize = cfg.IFREQ.winsize * header.samplingrate;
if winsize > size(data,2)
    winsize = size(data,2);
end
j = 1;
for i = Nact;
    if size(data,2) < 12000
        [SpectAllM(j,:),~,SpectAllF] = pmtm(data(i,:),4,winsize,header.samplingrate); %#ok<AGROW>
    else
        window = ones(1,winsize);
        [SpectAllM(j,:),SpectAllF] = pwelch(data(i,:),window,1,winsize,header.samplingrate); %#ok<AGROW>
    end
    j = j + 1;
end

freqlow = find(SpectAllF >= cfg.IFREQ.lowfreq, 1 );
freqhigh = find(SpectAllF <= cfg.IFREQ.highfreq, 1,'last');
deltafreq = round(1 / (SpectAllF(2)-SpectAllF(1)));
if ~isempty(freqlow) & ~isempty(freqhigh)
    Spect2 = SpectAllM(:,freqlow:freqhigh);
    [~,cc] = ndgrid(1:size(Spect2,1),1:size(Spect2,2));
    Spectt = sum(Spect2(:));
    c2 = sum(Spect2(:) .* cc(:)) / Spectt;
    header.IFREQ.MF = SpectAllF(freqlow -1) + (c2 * (SpectAllF(2) - SpectAllF(1)));
    iMF = round(((header.IFREQ.MF - SpectAllF(1))/(SpectAllF(2)-SpectAllF(1)))+1);
    [~,iPF] = max(SpectAllM(:,iMF-deltafreq:iMF+deltafreq),[],2);
    iPF = median(iPF);
    if iPF > 1 & iPF < (2*deltafreq+1)
        iPF = iPF + iMF - deltafreq - 1;
    else
        iPF = iMF;
    end
    header.IFREQ.PF = SpectAllF(iPF);
    freqlow2 = find(SpectAllF >= cfg.IFREQ.lowfreq-1, 1 );
    [~,iTF] = min(SpectAllM(:,freqlow2:iPF),[],2);
    iTF = iTF + freqlow2 - 1;
    if max(iTF) > 1
        iTF = round(median(iTF(iTF>1)));
        header.IFREQ.TF = SpectAllF(iTF);
    else
        header.IFREQ.TF = round(10*(SpectAllF(iMF) - SpectAllF(1) / 2)) / 10;
    end
else
    header.IFREQ = [];
end
AlphaCut = (header.IFREQ.TF + header.IFREQ.PF) / 2;
header.IFREQ.Bands = {'Delta',header.IFREQ.TF-4,header.IFREQ.TF-2,2,4; ...
                      'Theta',header.IFREQ.TF-2,header.IFREQ.TF,4,6; ...
                      'Alpha1',header.IFREQ.TF,header.IFREQ.PF,6,10; ...
                      'Alpha1a',header.IFREQ.TF,AlphaCut,6,8; ...
                      'Alpha1b',AlphaCut,header.IFREQ.PF,8,10; ...
                      'Alpha2',header.IFREQ.PF,header.IFREQ.PF+2,10,12; ...
                      'Beta',header.IFREQ.PF+2,header.IFREQ.PF+20,12,30; ...
                      'Beta1',header.IFREQ.PF+2,header.IFREQ.PF+11,12,21; ...
                      'Beta2',header.IFREQ.PF+11,header.IFREQ.PF+20,21,30; ...
                      'PF',header.IFREQ.PF-1,header.IFREQ.PF+1,9,11};

if ~isempty(IFREQ_file)
    IFREQ = header.IFREQ; %#ok<NASGU>
    save(IFREQ_file,'IFREQ');
end
if ~isempty(IFREQ_file)
    lab_plot_freqs(SpectAllM,SpectAllF,header.IFREQ,[IFREQ_file(1:end-6) '_IFREQ.jpg']);
    xlsout = cat(1,{'FreqBand','LowFreq','HighFreq'},header.IFREQ.Bands(:,1:3));
    lab_write_xls([IFREQ_file(1:end-6) '_IFREQ.xls'],xlsout);
end

return