% Spectral analysis
%           pwelch for data > 3 x FFT-window
%           pmtm   for data < 3 x FFT-window
%
% [Result,cfg] = lab_calculate_spectras(data,header,cfg)
%
% data = datamatrix (nchans x timeframes)
% header = Output of 'lab_read_data'
% cfg.Output_filepath
% cfg.Output_file
% cfg.EEG_filepath (optional)
%
% Written by F. Hatz 2012 Neurology Basel

function [Result,cfg] = lab_calculate_spectras(data,header,cfg)

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

if ~isfield(header,'goodchans')
    header.goodchans = 1:header.numdatachannels;
    header.badchans = [];
end

if isfield(cfg,'SKIP') & isfield(cfg.SKIP,'FFT') & cfg.SKIP.FFT == true;
    Result = [];
    return
end

disp('Spectral analysis')
if ~isfield(cfg,'FFT')
    [cfg,skipprocessing] = lab_set_calculate_spectras(cfg,header);
    if skipprocessing == 1
        Result = [];
        return
    else
        pause(0.1);
    end
end

winsize = cfg.FFT.winsize*header.samplingrate;
window = ones(1,winsize);
if cfg.FFT.hanpercent > 0
    hanpercent= cfg.FFT.hanpercent / 100;
    hanning = hann(ceil(winsize*hanpercent));
    overlap = hanpercent/2;
    window(1:ceil(overlap*winsize)) = hanning(1:ceil(overlap*winsize));
    window(end-ceil(overlap*winsize)+1:end) = hanning(end-ceil(overlap*winsize)+1:end);
else
    overlap = 0;
end
window = window';
clearvars hanning

if winsize > size(data,2)
    Result = [];
    disp('Abort FFT: segment to short for FFT-analysis')
    return
end

% Calculate new reference
if strcmp(cfg.FFT.eegsource,'montage') & isfield(cfg.FFT,'montage') & ~isempty(cfg.FFT.montage)
    [data,header,cfg.FFT] = lab_references(data,header,cfg.FFT.montage,cfg.FFT);
elseif ~strcmp(cfg.FFT.eegsource,'input')
    [data,header,cfg.FFT] = lab_references(data,header,cfg.FFT.eegsource,cfg.FFT);
end
Nact = header.numdatachannels;

% Reduce to data channels
if ~isfield(header,'numdatachannels')
    header.numdatachannels = size(data,1);
end
[data,header] = lab_reduce_channels(data,header,1:header.numdatachannels);

% interpolate bad
if isfield(cfg.FFT,'interpolatemethod') & ~isempty(cfg.FFT.interpolatemethod) & ~strcmp(cfg.FFT.interpolatemethod,'disabled')
    disp('   Interpolation of bad channels')
    [data,header] = lab_interpolate_bad(data,header,cfg.FFT.interpolatemethod,'noverbose');
end

% Define number of spectras
steps = 1:(winsize-fix(ceil(overlap*winsize))):(size(data,2)-winsize+1);

% look for bad channels
if length(steps) > 1 & isfield(cfg,'BADELEC') & isfield(cfg.BADELEC,'length') & ...
        cfg.BADELEC.length == (cfg.FFT.winsize - overlap*cfg.FFT.winsize) & ...
        isfield(header,'bad') & isfield(header.bad,'epochs') & ...
        size(header.bad.epochs,2) == size(steps,2)
    SpectValid = abs(header.bad.epochs - 1);
elseif length(steps) > 1
    disp('   Search for bad channels in epochs')
    settings = cfg.FFT.BAD;
    SpectValid = ones(size(data,1),size(steps,2));
    for step = 1:size(steps,2);
        datatmp = data(:,steps(1,step):steps(1,step)+winsize-1);
        settings.length = size(datatmp,2);
        bad = lab_detect_bad(datatmp,header,settings,cfg,'noverbose','silent');
        SpectValid(bad,step) = 0;
        clearvars bad
    end
else
    SpectValid = ones(size(data,1),1);
end

% Calculate Spectras
fprintf(['   Calculate spectras (method: ' cfg.FFT.method ')'])
stepcounter = floor(size(steps,2)/20);
if stepcounter < 1
    stepcounter = 1;
end
Cpower = [];
if strcmp(cfg.FFT.method,'wavelet') | strcmp(cfg.FFT.method,'rihaczek')
    if size(steps,2) > 1
        disp('    Take only first epoch (multiple epochs not supported for wavelet or rihaczek')
        datatmp = data(:,1:winsize);
    else
        datatmp = data;
    end
    if strcmp(cfg.FFT.method,'wavelet')
        % Define frequencies
        SpectF = cfg.FFT.lowfreq:header.samplingrate/winsize:cfg.FFT.highfreq;
        time = 0:1/header.samplingrate:winsize/header.samplingrate;
        time = time(1,1:winsize);
        Spect = lab_specest_wavelet(datatmp,time,'freqoi',SpectF);
        Spect = abs(permute(Spect,[3 2 1]));
    elseif strcmp(cfg.FFT.method,'rihaczek')
        deltafreq = header.samplingrate / winsize;
        lowfreq = round(cfg.FFT.lowfreq / deltafreq)+1;
        highfreq = round(cfg.FFT.highfreq / deltafreq)+1;
        SpectF = 0:deltafreq:header.samplingrate/2;
        SpectF = SpectF(1,lowfreq:highfreq);
        for i = 1:Nact;
            tmp = lab_rihaczek(datatmp(i,:)',header.samplingrate,winsize);
            if ~exist('Spect','var');
                Spect = abs(tmp(lowfreq:highfreq,:))';
            else
                Spect(:,:,i) = abs(tmp(lowfreq:highfreq,:))';
            end
        end
    end
    if isfield(cfg.FFT,'calcpower') & cfg.FFT.calcpower == true
        Cpower = datatmp.^2;
    end
    dotimefreq = true;
    disp(':')
else
    for step = 1:size(steps,2);
        datatmp = data(:,steps(1,step):steps(1,step)+winsize-1);
        for i = 1:Nact;
            if strcmp(cfg.FFT.method,'multitaper')
                [spectra,~,freqs] = pmtm(detrend(datatmp(i,:) .* window'),4,winsize,header.samplingrate);
            else
                [spectra,freqs] = pwelch(detrend(datatmp(i,:)),window,1,winsize,header.samplingrate);
            end
            if ~exist('Spect','var');
                deltafreq = freqs(2,1)-freqs(1,1);
                lowfreq = round(cfg.FFT.lowfreq / deltafreq)+1;
                highfreq = round(cfg.FFT.highfreq / deltafreq)+1;
                Spect = zeros(size(SpectValid,2),highfreq - lowfreq + 1,Nact);
                SpectF = freqs(lowfreq:highfreq,1)';
            end
            Spect(step,:,i) = spectra(lowfreq:highfreq)';
            clearvars spectra freqs;
        end
        if isfield(cfg.FFT,'calcpower') & cfg.FFT.calcpower == true
            Cpower(:,step) = sum(datatmp.^2,2) / size(datatmp,2); %#ok<AGROW>
        end
        clearvars datatmp headertmp
        if mod(step,stepcounter) == 0
            fprintf('.')
        end
    end
    disp(':')
    dotimefreq = false;
end

if size(Spect,1) > 1 & dotimefreq == false
    % Delete epochs with to mutch bad channels
    validtmp = [];
    if ~isfield(cfg.FFT,'percentgood')
        cfg.FFT.percentgood = 100;
        percentgood = 1;
    else
        percentgood = cfg.FFT.percentgood / 100;
    end
    while size(validtmp,2) < 1 | size(validtmp,2) < size(steps,2)/2
        validtmp = find(mean(SpectValid) >= percentgood);
        percentgood = percentgood - 0.01;
    end
    if cfg.FFT.percentgood ~= percentgood
        disp(['   Warning: not enough valid periods, PERCENT GOOD CHANNELS set to ' num2str((percentgood + 0.01)*100)])
    end
    
    % Delete non-valid ffts
    validtmp2 = find(max(max(isnan(Spect),[],3),[],2));
    validtmp = setdiff(validtmp,validtmp2);
    if isempty(validtmp)
        disp('   no valid epochs for spectra calculation found')
        Result = [];
        return
    end
    
    Spect = Spect(validtmp,:,:);
    if ~isempty(Cpower) & length(validtmp) == size(Cpower,2)
        Cpower = sqrt(mean(Cpower(:,validtmp),2));
    else
        Cpower = [];
    end
end

% Calculate mean, median, SD
Specttmp = permute(Spect,[3 2 1]).^0.5;
SpectMean = mean(Specttmp,3).^2; %#ok<NASGU>
SpectSDlow = (mean(Specttmp,3) - std(Specttmp,[],3)).^2; %#ok<NASGU>
SpectSDhigh = (mean(Specttmp,3) + std(Specttmp,[],3)).^2; %#ok<NASGU>
SpectMedian = median(Specttmp,3).^2; %#ok<NASGU>
if size(Spect,1) > 1 & dotimefreq == false
    SpectIncluded = zeros(1,size(SpectValid,2));
    SpectIncluded(validtmp) = 1;
else
    SpectIncluded = 1;
end
clearvars Specttmp validtmp

% Create output folder
if ~isfield(cfg.FFT,'folder')
    cfg.FFT.folder = 'Spect';
end
warning off %#ok<WNOFF>
mkdir (fullfile(cfg.Output_filepath,cfg.FFT.folder));
warning on %#ok<WNON>
Spect_filepath = fullfile(cfg.Output_filepath,cfg.FFT.folder);

if size(Spect,1) > 1 & dotimefreq == false
    %--------------------------------------------------------------------------
    % Write Marker Spectra file (*_Spectra.mrk)
    %--------------------------------------------------------------------------
    disp('   Write *_Spectra.mrk')
    SpectraMarker_file=fullfile(Spect_filepath,[cfg.Output_fileS '_Spectra.mrk']);
    Markers.POS = 1:(winsize-fix(ceil(overlap*winsize))):(size(data,2)-winsize+1);
    Markers.DUR = repmat(winsize,1,size(Markers.POS,2));
    Markers.TYP = cellstr(num2str((1:size(Markers.POS,2))'))';
    lab_write_mrk(SpectraMarker_file,Markers);
    clearvars Markers
end

%--------------------------------------------------------------------------
% Write Spectra summary file (Spectra.mat)
%--------------------------------------------------------------------------
cd(Spect_filepath);
disp('   Write Spectra.mat')
if dotimefreq == true & isfield(cfg,'EEG_filepath')
    [~,~,~,EEG_fileS] = lab_filename(cfg.EEG_file);
    SpectraAll_file = [EEG_fileS '_Spectra.mat'];
elseif isfield(cfg,'patient')
    SpectraAll_file = [cfg.patient '_Spectra.mat'];
elseif isfield(cfg,'EEG_filepath')
    tmp = strfind(cfg.EEG_filepath,filesep);
    SpectraAll_file = [cfg.EEG_filepath((tmp(end-1)+1):(end-1)) '_Spectra.mat'];
    clearvars tmp;
else
    SpectraAll_file = 'Spectra.mat';
end
if isfield(cfg.FFT,'deleteold') & islogical(cfg.FFT.deleteold) & cfg.FFT.deleteold == true
    warning off %#ok<WNOFF>
    if ~isfield(cfg,'listold')
        if exist(fullfile(Spect_filepath,SpectraAll_file),'file')
            disp('     delete Spectra.mat from previous run')
            delete(fullfile(Spect_filepath,SpectraAll_file));
        end
        cfg.listold = cellstr(fullfile(Spect_filepath,SpectraAll_file));
    elseif min(~strcmp(cfg.listold,fullfile(Spect_filepath,SpectraAll_file)))
        if exist(fullfile(Spect_filepath,SpectraAll_file),'file')
            disp('     delete Spectra.mat from previous run')
            delete(fullfile(Spect_filepath,SpectraAll_file));
        end
        cfg.listold = [cfg.listold cellstr(fullfile(Spect_filepath,SpectraAll_file))];
    end
    warning on %#ok<WNON>
end

% Combine results with previous ones
headertmp = header;
clearvars header
try %#ok<TRYNC>
    load(fullfile(Spect_filepath,SpectraAll_file));
end
if ~exist('SpectAll','var')
    SpectAll = [];
    SpectAllF = SpectF;
    SpectAllValid = [];
    SpectAllIncluded = {};
    header = headertmp;
else
    if isfield(header,'quality') & isfield(headertmp,'quality')
        header.quality = [header.quality headertmp.quality];
    end
    if isfield(header,'activationsexcluded') & isfield(headertmp,'activationsexcluded')
        header.activationsexcluded = [header.activationsexcluded headertmp.activationsexcluded];
    end
    if isfield(header,'badchans') & isfield(headertmp,'badchans')
        header.badchans = union(header.badchans,headertmp.badchans);
        if size(header.badchans,1) > 1
            header.badchans = header.badchans';
        end
        header.goodchans = setdiff(1:size(SpectAll,3),header.badchans); %#ok<NODEF>
        header.goodchans = header.goodchans(:)';
    end
    if isfield(header,'interpolated') & isfield(headertmp,'interpolated')
        header.interpolated = union(header.interpolated,headertmp.interpolated);
        if size(header.interpolated,1) > 1
            header.interpolated = header.interpolated';
        end
    end
end
SpectAll = cat(1,SpectAll,Spect);
SpectAllValid = [SpectAllValid SpectValid]; %#ok<NASGU>
SpectAllIncluded{1,end+1} = SpectIncluded; %#ok<NASGU>
if size(SpectAll,1) > 1
    for i = 1:Nact;
        SpectAllMean(i,:) = mean(SpectAll(:,:,i).^0.5).^2; %#ok<AGROW>
        SpectAllSDlow(i,:)  = (mean(SpectAll(:,:,i).^0.5) - std(SpectAll(:,:,i).^0.5)).^2; %#ok<AGROW,NASGU>
        SpectAllSDhigh(i,:)  = (mean(SpectAll(:,:,i).^0.5) + std(SpectAll(:,:,i).^0.5)).^2; %#ok<AGROW,NASGU>
        SpectAllMedian(i,:)  = median(SpectAll(:,:,i)); %#ok<AGROW>
    end
else
    SpectAllMean = permute(SpectAll,[3 2 1]);
    SpectAllSDlow  = zeros(Nact,size(SpectAll,2)); %#ok<NASGU>
    SpectAllSDhigh  = zeros(Nact,size(SpectAll,2)); %#ok<NASGU>
    SpectAllMedian = permute(SpectAll,[3 2 1]);
end
if ~exist('CpowerAll','var')
    CpowerAll = [];
end
if ~isempty(Cpower)
    CpowerAll = [CpowerAll Cpower]; %#ok<NASGU>
end
save(fullfile(Spect_filepath,SpectraAll_file),'SpectAll','SpectAllF','SpectAllValid', ...
    'SpectAllMean','SpectAllSDlow','SpectAllSDhigh','SpectAllMedian','SpectAllIncluded','CpowerAll','header','dotimefreq');

% Store filepath to Spectra.mat in cfg
if ~isfield(cfg.FFT,'filepaths')
    cfg.FFT.filepaths{1,1} = fullfile(Spect_filepath,SpectraAll_file);
elseif max(strcmp(cfg.FFT.filepaths,fullfile(Spect_filepath,SpectraAll_file))) == 0
    cfg.FFT.filepaths = [cfg.FFT.filepaths cellstr(fullfile(Spect_filepath,SpectraAll_file))];
end

% Collect results
Result.Power_mean = SpectAllMean;
Result.Power_median = SpectAllMedian;
Result.Power_freqs = SpectAllF;

% Save xls files
if (~isfield(cfg,'lastsegment') | cfg.lastsegment == true) & ...
        (~isfield(cfg.FFT,'writexls') | cfg.FFT.writexls == true)
    cd(Spect_filepath);
    ChannelsTitel = [{'Channel'} cellstr(header.channels(1:Nact,:))'];
    Spectrum = num2cell(SpectAllMean);
    SpectrumF = num2cell(SpectAllF);
    Spectrum = [SpectrumF' Spectrum']';
    Spectrum = cat(2,ChannelsTitel',Spectrum);
    if size(Spectrum,2) > 255
        fileout = [SpectraAll_file(1:end-4) 'Mean.xlsx'];
    else
        fileout = [SpectraAll_file(1:end-4) 'Mean.xls'];
    end
    lab_write_xls(fileout, Spectrum);
    Spectrum = num2cell(SpectAllMedian);
    Spectrum = [SpectrumF' Spectrum']';
    Spectrum = cat(2,ChannelsTitel',Spectrum);
    if size(Spectrum,2) > 255
        fileout = [SpectraAll_file(1:end-4) 'Median.xlsx'];
    else
        fileout = [SpectraAll_file(1:end-4) 'Median.xls'];
    end
    lab_write_xls(fileout,Spectrum);
    
    fileout = [SpectraAll_file(1:end-4) 'ChannelPower.xlsx'];
    xlsout = [cellstr(header.channels(1:Nact,:)) num2cell(sqrt(mean(Cpower.^2,2)))];
    lab_write_xls(fileout,xlsout);
    
    %--------------------------------------------------------------------------
    % Calculate band power
    %--------------------------------------------------------------------------
    disp('   Calculate BandPower')
    if isfield(cfg.FFT,'spectralbandsI') & cfg.FFT.spectralbandsI == true
        if ~isfield(header,'IFREQ') | ~isfield(header.IFREQ,'Bands') | isempty(header.IFREQ.Bands)
            [~,header] = lab_indiv_freqbands(data,header);
        end
        if ~isempty(header.IFREQ.Bands)
            spectralbands = cell2mat(header.IFREQ.Bands(:,2:5));
        else
            disp('   Abort: Calculation of individual bands not possible')
            return
        end
    else
        spectralbands = cell2mat(cfg.FFT.spectralbands(:,2:3));
        spectralbands = [spectralbands spectralbands];
    end
    BandTitel = cell(1,size(spectralbands,1));
    for i = 1:size(spectralbands,1);
        BandTitel{1,i} = ['F' num2str(spectralbands(i,3)) 'F' num2str(spectralbands(i,4))];
    end
    
    %mean
    BandPower = zeros(Nact,size(spectralbands,1));
    for i = 1:size(spectralbands,1);
        range = spectralbands(i,1):.01:spectralbands(i,2);
        for j = 1:size(SpectAllMean,1)
            int_spec = exp(interp1(log(SpectAllF),log(SpectAllMean(j,:)),log(range),'linear','extrap'));
            BandPower(j,i) = trapz(range,int_spec);
        end
    end
    Result.BandPower_mean = BandPower;
    Result.BandPower_freqs = BandTitel;
    BandPower = num2cell(BandPower);
    BandPower = [BandTitel' BandPower']';
    BandPower = cat(2,ChannelsTitel',BandPower);
    if size(BandPower,2) > 255
        fileout = [SpectraAll_file(1:end-11) 'BandPowerMean.xlsx'];
    else
        fileout = [SpectraAll_file(1:end-11) 'BandPowerMean.xls'];
    end
    lab_write_xls(fileout, BandPower);
    clearvars range int_spec i j BandPower BandTitel
    
    %median
    BandPower = zeros(Nact,size(spectralbands,1));
    for i = 1:size(spectralbands,1);
        range = spectralbands(i,1):.01:spectralbands(i,2);
        for j = 1:size(SpectAllMedian,1)
            int_spec = exp(interp1(log(SpectAllF),log(SpectAllMedian(j,:)),log(range),'linear','extrap'));
            BandPower(j,i) = trapz(range,int_spec);
        end
    end
    Result.BandPower_median = BandPower;
    BandPower=num2cell(BandPower);
    BandPower = [BandTitel' BandPower']';
    BandPower = cat(2,ChannelsTitel',BandPower);
    if size(BandPower,2) > 255
        fileout = [SpectraAll_file(1:end-11) 'BandPowerMedian.xlsx'];
    else
        fileout = [SpectraAll_file(1:end-11) 'BandPowerMedian.xls'];
    end
    lab_write_xls(fileout, BandPower);
    clearvars range int_spec i j BandPower BandTitel
    
    if dotimefreq == true
        OutputFreqs_filepath = fullfile(Spect_filepath,'TimeFreqs');
        warning off %#ok<WNOFF>
        mkdir (OutputFreqs_filepath);
        warning on %#ok<WNON>
        f = figure('Visible','off');
        for i = 1:Nact
            imagesc(Spect(:,:,i)');
            set(gca,'YTickLabel',num2str(SpectF(get(gca,'YTick'))'),'FontName','Times','fontsize',9);
            title([SpectraAll_file(1:end-11) '_TimeFreqs' num2str(i) '.jpg']);
            lab_print_figure(fullfile(OutputFreqs_filepath,[SpectraAll_file(1:end-11) '_TimeFreqs' num2str(i) '.jpg']),f);
            clf;
            dlmwrite(fullfile(OutputFreqs_filepath,[SpectraAll_file(1:end-11) '_TimeFreqs' num2str(i) '.csv']),Spect(:,:,i)');
        end
        close(f);
        dlmwrite(fullfile(OutputFreqs_filepath,[SpectraAll_file(1:end-11) '_Freqs.csv']),SpectF','delimiter',';');
    end
    
    %--------------------------------------------------------------------------
    % Calculate band power for Mappings
    %--------------------------------------------------------------------------
    if isfield(cfg.FFT,'Mappings') & ~isempty(cfg.FFT.Mappings)
        if cfg.FFT.Mappings.mappingsChannels ~= size(data,1)
            cfg.FFT.Mappings = lab_reduce_mappings(cfg.FFT.Mappings,[],cfg);
        end
        if cfg.FFT.Mappings.mappingsChannels == size(data,1)
            disp('   Calculate BandPower for Mappings')
            Mappings = cfg.FFT.Mappings;
            if strcmp(Mappings.mappingstitle(end,1),'no region')
                Mappings.mappings = Mappings.mappings(1,1:end-1);
                Mappings.mappingstitle = Mappings.mappingstitle(1:end-1,1);
            end
                        
            % Calculate Mappings Frequencies and plot/save Spectras
            for i = 1 : size(Mappings.mappings,2);
                SpectMappingsMean(i,:)= mean(SpectAllMean(Mappings.mappings{i},:).^0.5).^2; %#ok<AGROW>
                SpectMappingsMedian(i,:)= median(SpectAllMedian(Mappings.mappings{i},:)); %#ok<AGROW>
            end
            
            if isfield(cfg.FFT,'plotspectras') & cfg.FFT.plotspectras == true
                OutputPlots_filepath = fullfile(Spect_filepath,'Plots');
                warning off %#ok<WNOFF>
                mkdir (OutputPlots_filepath);
                warning on %#ok<WNON>
                
                freqlabel = cell(1,length(SpectAllF));
                for i = 1:length(SpectAllF)
                    if round(SpectAllF(i)) == SpectAllF(i)
                        freqlabel(i) = num2cell(SpectAllF(i));
                    end
                end
                tmp = find(~cellfun(@isempty,freqlabel));
                if ~isempty(tmp) & length(tmp) >25
                    for i = 2:2:length(tmp)
                        freqlabel{1,tmp(i)} = [];
                    end
                end
                
                f = figure('Visible','off');
                for i = 1 : size(Mappings.mappings,2);
                    area(SpectMappingsMean(i,:));
                    set(gca,'XTick',1:length(freqlabel));
                    set(gca,'XTickLabel',freqlabel,'FontName','Times','fontsize',9);
                    title(Mappings.mappingstitle(i,1));
                    lab_print_figure(fullfile(OutputPlots_filepath,[SpectraAll_file(1:end-11) '_' Mappings.mappingstitle{i} '_mean.jpg']),f);
                    clf;
                    
                    area(SpectMappingsMedian(i,:));
                    set(gca,'XTick',1:length(freqlabel));
                    set(gca,'XTickLabel',freqlabel,'FontName','Times','fontsize',9);
                    title(Mappings.mappingstitle(i,1));
                    lab_print_figure(fullfile(OutputPlots_filepath,[SpectraAll_file(1:end-11) '_' Mappings.mappingstitle{i} '_median.jpg']),f);
                    clf;
                end
                close(f);
            end
            
            %mean
            RegionTitel = [{'Region'} Mappings.mappingstitle'];
            BandPower = zeros(size(Mappings.mappings,2),size(spectralbands,1));
            for i = 1:size(spectralbands,1);
                range = spectralbands(i,1):.01:spectralbands(i,2);
                for j = 1:size(SpectMappingsMean,1)
                    int_spec = exp(interp1(log(SpectAllF),log(SpectMappingsMean(j,:)),log(range),'linear','extrap'));
                    BandPower(j,i) = trapz(range,int_spec);
                end
            end
            Result.BP_Mapping_mean = BandPower;
            BandPower=num2cell(BandPower);
            BandPower = [BandTitel' BandPower']';
            BandPower = cat(2,RegionTitel',BandPower);
            if size(BandPower,2) > 255
                fileout = [SpectraAll_file(1:end-11) 'BandPowerMappingMean.xlsx'];
            else
                fileout = [SpectraAll_file(1:end-11) 'BandPowerMappingMean.xls'];
            end
            lab_write_xls(fileout, BandPower);
            clearvars range int_spec i j BandPower BandTitel
            
            %median
            BandPower = zeros(size(Mappings.mappings,2),size(spectralbands,1));
            for i = 1:size(spectralbands,1);
                range = spectralbands(i,1):.01:spectralbands(i,2);
                for j = 1:size(SpectMappingsMedian,1)
                    int_spec = exp(interp1(log(SpectAllF),log(SpectMappingsMedian(j,:)),log(range),'linear','extrap'));
                    BandPower(j,i) = trapz(range,int_spec);
                end
            end
            Result.BP_Mapping_median = BandPower;
            BandPower=num2cell(BandPower);
            BandPower = [BandTitel' BandPower']';
            BandPower = cat(2,RegionTitel',BandPower);
            if size(BandPower,2) > 255
                fileout = [SpectraAll_file(1:end-11) 'BandPowerMappingMedian.xlsx'];
            else
                fileout = [SpectraAll_file(1:end-11) 'BandPowerMappingMedian.xls'];
            end
            lab_write_xls(fileout, BandPower);
            clearvars range int_spec i j BandPower BandTitel
        else
            disp('    Mappings not calculated (mismatch mappings-file - data)')
        end
    end
end

%--------------------------------------------------------------------------
% Write verbose file (*.vrb)
%--------------------------------------------------------------------------
if strcmp(SpectraAll_file,'Spectra.mat')
    fid=fopen(fullfile(Spect_filepath,'Spectra.vrb'),'w');
else
    fid=fopen(fullfile(Spect_filepath,[SpectraAll_file(1:end-12) '_Spectra.vrb']),'w');
end
fprintf(fid,'Power analysis\n');
fprintf(fid,['Patient: ' SpectraAll_file(1:end-11)]);
fprintf(fid,'\n');
fprintf(fid,['EEG reference: ' cfg.FFT.eegsource]);
fprintf(fid,'\n');
fprintf(fid,['FFT Window: ' num2str(winsize) ' with ' num2str(cfg.FFT.hanpercent) '% Hanning-Window']);
fprintf(fid,'\n');
fprintf(fid,['Filter: ' num2str(roundn(SpectF(1,1),2)) 'Hz  - ' num2str(roundn(SpectF(1,end),2)) 'Hz']);
fprintf(fid,'\n');
fprintf(fid,['Method: ' cfg.FFT.method]);
fprintf(fid,'\n');
if isfield(cfg.FFT,'interpolatemethod') & ~strcmp(cfg.FFT.interpolatemethod,'disabled')
    fprintf(fid,['Bad electrodes interpolated: ' cfg.FFT.interpolatemethod]);
end
fprintf(fid,'\n');
if isfield(cfg.FFT,'BAD')
    fprintf(fid,'Settings for detecting bad channels:\n');
    fprintf(fid,['  Bad Spectra 50Hz (Limit ' num2str(round(cfg.FFT.BAD.freqlim50)) 'percent)\n']);
    fprintf(fid,['  Bad Spectra 60Hz (Limit ' num2str(round(cfg.FFT.BAD.freqlim60)) 'percent)\n']);
    fprintf(fid,['  Bad Spectra '  num2str(cfg.FFT.BAD.spectshigh(1,1)) '-' num2str(cfg.FFT.BAD.spectshigh(1,2)) ...
        'Hz (Limit ' num2str(round(cfg.FFT.BAD.freqlimhigh)) 'percent)\n']);
    fprintf(fid,['  Bad Spectra '  num2str(cfg.FFT.BAD.spectslow(1,1)) '-' num2str(cfg.FFT.BAD.spectslow(1,2)) ...
        'Hz (Limit ' num2str(round(cfg.FFT.BAD.freqlimlow)) 'percent)\n']);
    fprintf(fid,['  Bad variance (zValue: '  num2str(cfg.FFT.BAD.zvaluevars) ')\n']);
    fprintf(fid,['  Bad hurst (zValue: '  num2str(cfg.FFT.BAD.zvaluehurst) ')\n']);
    if exist('percentgood','var')
        fprintf(fid,['  Percent of good channels per frequency window: '  num2str(percentgood*100) '\n']);
    end
end
if exist('BandTitel','var')
    fprintf(fid,'Spectral bands\n');
    tmp = char(BandTitel);
    tmp = [tmp char(32*ones(size(tmp,1), 3))]';
    fprintf(fid,tmp);
end
fclose(fid);

return