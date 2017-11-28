% Helper function for lab_collect_spectras
%
% [cfg,result,calc,Mappings,skipprocessing] = lab_collect_spectraldata_header(cfg,calc,result,Mappings)
%
% written by F. Hatz 2012

function [cfg,result,calc,Mappings,skipprocessing] = lab_collect_spectraldata_header(cfg,calc,result,Mappings)

skipprocessing = 0;

if isfield(calc,'SpectAll')
    SpectAll = calc.SpectAll;
    SpectAllF = calc.SpectAllF;
else
    load(calc.matfile);
end
if isfield(cfg.CollectFFT,'spectralbandsI') & cfg.CollectFFT.spectralbandsI == true & isfield(header,'IFREQ') & isfield(header.IFREQ,'Bands')
    Spectralbands = cell2mat(header.IFREQ.Bands(:,2:5));
    if isnumeric(cfg.CollectFFT.spectralbands) & ~isempty(cfg.CollectFFT.spectralbands)
        Spectralbands = Spectralbands(cfg.CollectFFT.spectralbands,:);
    end
elseif ~isempty(cfg.CollectFFT.spectralbands) & size(cfg.CollectFFT.spectralbands,2) >= 3
    Spectralbands = cell2mat(cfg.CollectFFT.spectralbands(:,2:3));
    Spectralbands = [Spectralbands Spectralbands];
else
    disp('Abort: Error with spectral bands')
    skipprocessing = 1;
    return
end
if exist('SpectAll','var') & exist('SpectAllF','var')
    cfg.numfreqbins = size(SpectAllF,2);
    if isfield(Mappings,'mappingsChannels') & size(SpectAll,3) ~= Mappings.mappingsChannels
        Mappings = lab_reduce_mappings(Mappings,cfg);
    end
    if isfield(Mappings,'mappingsChannels') & size(SpectAll,3) ~= Mappings.mappingsChannels
        Mappings = [];
        disp('Mappings disabled - channels not matching')
    end
    calc.numchannels = size(SpectAll,3);
    if isnumeric(SpectAllF)
        calc.bandtitles = cell(1,size(Spectralbands,1));
        for i = 1:(size(Spectralbands,1));
            n = find(SpectAllF >= Spectralbands(i,3), 1 );
            m = find(SpectAllF < Spectralbands(i,4), 1, 'last' );
            calc.bandtitles{1,i} = ['F' num2str(round(SpectAllF(1,n)*100)/100) ':' num2str(round(SpectAllF(1,m+1)*100)/100)];
        end
    else
        calc.bandtitles = SpectAllF;
    end
    
    % Calculate header Mapping
    if exist('Mappings','var') & ~isempty(Mappings)
        resultheader = cell(1,(size(calc.bandtitles,2)*size(Mappings.mappings,2)));
        for j = 1:size(calc.bandtitles,2)
            for i = 1:size(Mappings.mappings,2)
                if isfield(Mappings,'mappingstitleS')
                    resultheader{1,((j-1)*size(Mappings.mappings,2)+i)} = [calc.bandtitles{1,j} '_' Mappings.mappingstitleS{i}];
                else
                    resultheader{1,((j-1)*size(Mappings.mappings,2)+i)} = [calc.bandtitles{1,j} '_' 'M' num2str(i)];
                end
            end
        end
        resultheader = [cellstr('MapFreq') resultheader];
        result.BP = resultheader;
        result.BPR = resultheader;
        result.BPM = resultheader;
        result.BPMR = resultheader;
        result.BPlog = resultheader;
        result.BPRlog = resultheader;
        result.BPMlog = resultheader;
        result.BPMRlog = resultheader;
        result.BE = resultheader;
        result.BEM = resultheader;   
        clearvars resultheader j i m n
        
        resultheader = cell(1,size(Mappings.mappings,2));
        for i = 1:size(Mappings.mappings,2)
            resultheader{1,i} = ['Median_Freq_M' num2str(i)];
        end
        resultheader = resultheader(:)';
        resultheader = [cellstr('MapFreq') resultheader];
        result.MF = resultheader;
        
        resultheader = cell(1,size(Mappings.mappings,2));
        for i = 1:size(Mappings.mappings,2)
            resultheader{1,i} = ['Peak_Freq_M' num2str(i)];
        end
        resultheader = resultheader(:)';
        resultheader = [cellstr('MapFreq') resultheader];
        result.PF = resultheader;
        
        resultheader = cell(1,size(Mappings.mappings,2));
        for i = 1:size(Mappings.mappings,2)
            resultheader{1,i} = ['Entropy_M' num2str(i)];
        end
        resultheader = resultheader(:)';
        resultheader = [cellstr('MapFreq') resultheader];
        result.TE = resultheader;
        
        resultheader = cell(1,size(Mappings.mappings,2));
        for i = 1:size(Mappings.mappings,2)
            resultheader{1,i} = ['EntropyM_M' num2str(i)];
        end
        resultheader = resultheader(:)';
        resultheader = [cellstr('MapFreq') resultheader];
        result.TEM = resultheader;
    end
    
    % Calculate header Every Channel
    resultheaderE = cell(1,(size(calc.bandtitles,2)*size(SpectAll,3)));
    for j = 1:size(calc.bandtitles,2)
        for i = 1:size(SpectAll,3)
            if exist('header','var') & isfield(header,'channels') & size(SpectAll,3) == size(header.channels,1)
                tmp = textscan(header.channels(i,:),'%s');
                tmp = tmp{1,1}{1,1};
                resultheaderE{1,((j-1)*size(SpectAll,3)+i)} = [calc.bandtitles{1,j} '_' tmp];
                clearvars tmp
            else
                resultheaderE{1,((j-1)*size(SpectAll,3)+i)} = [calc.bandtitles{1,j} '_' 'E' num2str(i)];
            end
        end
    end
    resultheaderE = [cellstr('MapFreq') resultheaderE];
    result.EBP = resultheaderE;
    result.EBPR = resultheaderE;
    result.EBPM = resultheaderE;
    result.EBPMR = resultheaderE;
    result.EBPlog = resultheaderE;
    result.EBPRlog = resultheaderE;
    result.EBPMlog = resultheaderE;
    result.EBPMRlog = resultheaderE;
    result.EBE = resultheaderE;
    result.EBEM = resultheaderE;
    clearvars resultheaderE i j
    
    % Calculate header Peak and Cog Freq
    result.BA = [cellstr(' ') cellstr(['Median_Freq_' num2str(cfg.CollectFFT.lowfreqcog) '_' num2str(cfg.CollectFFT.highfreqcog) 'Hz']) ...
        cellstr(['Peak_Freq_' num2str(cfg.CollectFFT.lowfreqpeak) '_' num2str(cfg.CollectFFT.highfreqpeak) 'Hz'])];
    result.Q = [cellstr(' ') cellstr(['Median_Freq_' num2str(cfg.CollectFFT.lowfreqcog) '_' num2str(cfg.CollectFFT.highfreqcog) 'Hz']) ...
        cellstr(['BandPower_' num2str(cfg.CollectFFT.lowfreqcog) '_' num2str(cfg.CollectFFT.highfreqcog) 'Hz']) ...
        cellstr(['Peak_Freq_' num2str(cfg.CollectFFT.lowfreqpeak) '_' num2str(cfg.CollectFFT.highfreqpeak) 'Hz']) ...
        cellstr(['BandPower_' num2str(cfg.CollectFFT.lowfreqpeak) '_' num2str(cfg.CollectFFT.highfreqpeak) 'Hz']) ...
        cellstr('Peak_Amplitude') cellstr('Peak2Min') cellstr('Peak_Ratio')];
    if ~isfield(cfg.CollectFFT,'BandTitel')
        result.Q = [result.Q  cellstr(['BandPower_' num2str(Spectralbands(end,3)) '_' num2str(Spectralbands(end,4)) 'Hz'])];
    else
        result.Q = [result.Q  cellstr('BandPower_i2_i30Hz')];
    end
    numchans = 1:size(SpectAll,3);
    if exist('header','var') & isfield(header,'channels') & size(SpectAll,3) == size(header.channels,1)
        result.EPF = [cellstr(' ') cellstr(header.channels(1:size(SpectAll,3),:))'];
    else
        result.EPF = [cellstr(' ') num2cell(numchans)];
    end
    result.EMF = result.EPF;
    result.ETE = result.EPF;
    result.ETEM = result.EPF;
    clearvars resultheader resultheaderE m n allchans exclude;
    
    
    if ~isfield(calc,'freqlabel')
        if iscell(SpectAllF)
            calc.freqlabel = SpectAllF;
            calc.freqlabelall = SpectAllF;
        else
            calc.freqlabel = cell(1,length(SpectAllF));
            calc.freqlabelall = num2cell(SpectAllF);
            for i = 1:length(SpectAllF)
                if round(SpectAllF(i)) == SpectAllF(i)
                    calc.freqlabel(i) = num2cell(SpectAllF(i));
                end
            end
            tmp = find(~cellfun(@isempty,calc.freqlabel));
            if ~isempty(tmp) & length(tmp) >25
                for i = 2:2:length(tmp)
                    calc.freqlabel{1,tmp(i)} = [];
                end
            elseif isempty(tmp) | length(tmp) < 5
                calc.freqlabel = cell(1,length(SpectAllF));
                for i = 1:floor(length(SpectAllF)/25):length(SpectAllF)
                    calc.freqlabel{1,i} = SpectAllF(i);
                end
            end
            clearvars tmp
        end
    end
else
    skipprocessing = 1;
    disp('Abort - calculation of header information not possible, wrong input data')
end