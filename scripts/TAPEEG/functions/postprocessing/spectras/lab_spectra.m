% Spectral analysis
%           pwelch for data > 3 x FFT-window
%           pmtm   for data < 3 x FFT-window
%
% [Result,cfg] = lab_spectra(data,header,cfg)
%
% data = datamatrix (nchans x timeframes)
% header = Output of 'lab_read_data'
% cfg.Output_filepath
% cfg.Output_file
% cfg.EEG_filepath (optional)
%
% Written by F. Hatz 2014 Neurology Basel

function [Result,cfg] = lab_spectra(data,header,cfg)

disp('Spectral analysis')
skipprocessing = 0;

Result = [];

if ~exist('cfg','var')
    cfg = [];
end
if ~isfield(cfg,'FFT')
    [cfg,skipprocessing] = lab_set_spectra(cfg);
    pause(0.1);
    if skipprocessing == 1
        return
    end
end

if ~exist('header','var')
    header = lab_create_header(data);
end
if ~isfield(header,'goodchans')
    header.goodchans = 1:header.numdatachannels;
    header.badchans = [];
end

if skipprocessing == 0
    winsize=cfg.FFT.winsize*header.samplingrate;
    if winsize > size(data,2) | cfg.FFT.winsize == 0;
        winsize = size(data,2);
    end
    if isfield(cfg.FFT,'hanpercent')
        hanpercent = cfg.FFT.hanpercent / 100;
    else
        hanpercent = 0;
    end
    hanning=hann(ceil(winsize * hanpercent));
    overlap = hanpercent / 2;
    window=ones(1,winsize);
    window(1:ceil(overlap*winsize))=hanning(1:ceil(overlap*winsize));
    window(end-ceil(overlap*winsize)+1:end)=hanning(end-ceil(overlap*winsize)+1:end);
    window=window';
    Nact=header.numdatachannels;
    clearvars hanning

    % Reduce to data channels
    if ~isfield(header,'numdatachannels')
        header.numdatachannels = size(data,1);
    end
    [data,header] = lab_reduce_channels(data,header,1:header.numdatachannels);
    
    % Calculate Spectras
    fprintf(['   Calculate spectras (method: ' cfg.FFT.method ')'])
    %[spectra,spectraP,SpectF,spectraV] = lab_pwelch(data(i,:),window,overlap,length(window),header.samplingrate,0.95,cfg.FFT.SD,'long-mean','powerarray');
    %[spectra,~,SpectF] = lab_pwelch(data(i,:),window,overlap,length(window),header.samplingrate,0.95,'linear','powerarray');
    steps = 1:(winsize-fix(ceil(overlap*winsize))):(size(data,2)-winsize+1);
    if size(steps,2) == 1
        cfg.FFT.method = 'multitaper';
    end
    for step = 1:size(steps,2);
        datatmp = data(:,steps(1,step):steps(1,step)+winsize-1);
        for i = 1:Nact;
            if strcmp(cfg.FFT.method,'multitaper')
                [spectra,~,freqs] = pmtm(detrend(datatmp(i,:)),4,winsize,header.samplingrate);
            else
                [spectra,freqs] = pwelch(detrend(datatmp(i,:)),window,1,winsize,header.samplingrate);
            end
            if ~ exist('Spect','var');
                deltafreq = freqs(2,1)-freqs(1,1);
                lowfreq = round(cfg.FFT.lowfreq / deltafreq)+1;
                highfreq = round(cfg.FFT.highfreq / deltafreq)+1;
                Spect = zeros(Nact,highfreq - lowfreq + 1,size(steps,2));
                SpectF=freqs(lowfreq:highfreq,1)';
            end
            Spect(i,:,step)=spectra(lowfreq:highfreq)';
            clearvars spectra freqs;
        end
        clearvars datatmp headertmp
        fprintf('.')
    end
    disp(' ')
    
    % Delete non-valid ffts
    tmp = find(max(max(isnan(Spect),[],2),[],1));
    validtmp = setdiff(1:size(Spect,3),tmp);
    clearvars tmp
    if isempty(validtmp)
        disp('   no valid spectra found')
        skipprocessing = 1;
    end
end
if skipprocessing == 0
    Spect = Spect(:,:,validtmp);
    % Calculate mean, median, SD
    SpectMean = mean(mean(Spect.^0.5,3),1).^2;
    SpectSDlow = (mean(mean(Spect.^0.5,3),1) - mean(std(Spect.^0.5,[],3),1)).^2;
    SpectSDhigh = (mean(mean(Spect.^0.5,3),1) + mean(std(Spect.^0.5,[],3),1)).^2;
    SpectMedian = median(median(Spect,3),1);
    SpectMedianLow = zeros(size(SpectMedian));
    SpectMedianHigh = zeros(size(SpectMedian));
    tmp = median(Spect,3);
    for i = 1:size(Spect,2)
        SpectMedianLow(1,i) = median(tmp(tmp(:,i) <= SpectMedian(i),i));
        SpectMedianHigh(1,i) = median(tmp(tmp(:,i) >= SpectMedian(i),i)); 
    end
    clearvars tmp
    clearvars validtmp
    Result.SpectMean = SpectMean;
    Result.SpectSDlow = SpectSDlow;
    Result.SpectSDhigh = SpectSDhigh;
    Result.SpectMedian = SpectMedian;
    Result.SpectMedianLow = SpectMedianLow;
    Result.SpectMedianHigh = SpectMedianHigh;
    Result.Spect = Spect;
    Result.Freqs = SpectF;
    
    freqlabel = cell(1,length(SpectF));
    for i = 1:length(SpectF)
        if round(SpectF(i)) == round(10*SpectF(i))/10 & mod(round(SpectF(i)),cfg.FFT.freqstep) == 0
            freqlabel(i) = num2cell(SpectF(i));
        end
    end
    if isfield(cfg.FFT,'average') & strcmp(cfg.FFT.average,'mean')
        figure('Color',[1 1 1]);
        h = area(SpectMean);
        set(h(1),'FaceColor',[0.6 0.6 0.6]);
        set(gca,'XTick',1:length(freqlabel),'XTickLabel',freqlabel,'FontName','Times','fontsize',9,'Box','off');
        title(['Spectral density plot ' num2str(cfg.FFT.lowfreq) '-' num2str(cfg.FFT.highfreq) 'Hz (mean +- std)']);
        hold on
        plot(SpectSDhigh,'--k');
        plot(SpectSDlow,'--k');
    else
        figure('Color',[1 1 1]);
        h = area(SpectMedian);
        set(h(1),'FaceColor',[0.6 0.6 0.6]);
        set(gca,'XTick',1:length(freqlabel),'XTickLabel',freqlabel,'FontName','Times','fontsize',9,'Box','off');
        title(['Spectral density plot ' num2str(cfg.FFT.lowfreq) '-' num2str(cfg.FFT.highfreq) 'Hz (median +- quartile)']);
        hold on
        plot(SpectMedianHigh,'--k');
        plot(SpectMedianLow,'--k');
    end
end

return