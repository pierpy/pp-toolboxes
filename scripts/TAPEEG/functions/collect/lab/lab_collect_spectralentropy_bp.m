% Helper function for lab_collect_spectras
%
% [cfg,result] = lab_collect_spectralentropy_bp(cfg,result,calc,SpectAll,SpectAllMean,SpectAllMedian,SpectAllF,Mappings)
%
% written by F. Hatz 2016

function [cfg,result] = lab_collect_spectralentropy_bp(cfg,result,calc,SpectAll,SpectAllMean,SpectAllMedian,SpectAllF,Mappings)

if isempty(SpectAllMean)
    emptymean = 1; %#ok<NASGU>
    SpectAllMean = permute(mean(SpectAll.^0.5,1).^2,[3 2 1]);
end
if isempty(SpectAllMedian)
    emptymedian = 1; %#ok<NASGU>
    SpectAllMedian = permute(median(SpectAll,1),[3 2 1]);
end

% ---------------------------------
% Calculate BandPower Every Channel
% ---------------------------------
entropyE = zeros(size(SpectAll,3),size(cfg.CollectFFT.freqs,1));
entropyEM = zeros(size(SpectAll,3),size(cfg.CollectFFT.freqs,1));
for i = 1:(size(cfg.CollectFFT.freqs,1));
    range = cfg.CollectFFT.freqs(i,1):.01:cfg.CollectFFT.freqs(i,2);
    for j = 1:size(SpectAllMean,1)
        int_spec = exp(interp1(log(SpectAllF),log(SpectAllMean(j,:)),log(range),'linear','extrap'));
        entropyE(j,i) = lab_spectral_entropy(int_spec);
        int_spec = exp(interp1(log(SpectAllF),log(SpectAllMedian(j,:)),log(range),'linear','extrap'));
        entropyEM(j,i) = lab_spectral_entropy(int_spec);
    end
end
clearvars j range int_spec

% calculate total power
range = min(cfg.CollectFFT.freqs(:,1)):.01:max(cfg.CollectFFT.freqs(:,2));
for j = 1:size(SpectAllMean,1)
    int_spec = exp(interp1(log(SpectAllF),log(SpectAllMean(j,:)),log(range),'linear','extrap'));
    entropyES(j,1) = lab_spectral_entropy(int_spec); %#ok<AGROW>
    int_spec = exp(interp1(log(SpectAllF),log(SpectAllMedian(j,:)),log(range),'linear','extrap'));
    entropyEMS(j,1) = lab_spectral_entropy(int_spec); %#ok<AGROW>
end
clearvars j range int_spec

% reshape data to output format
entropyE = reshape(entropyE,1,size(entropyE,1)*size(entropyE,2));
entropyES = reshape(entropyES,1,size(entropyES,1)*size(entropyES,2));
entropyE=[cellstr(calc.patient) num2cell(entropyE)];
entropyES=[cellstr(calc.patient) num2cell(entropyES)];
entropyEM = reshape(entropyEM,1,size(entropyEM,1)*size(entropyEM,2));
entropyEMS = reshape(entropyEMS,1,size(entropyEMS,1)*size(entropyEMS,2));
entropyEM=[cellstr(calc.patient) num2cell(entropyEM)];
entropyEMS=[cellstr(calc.patient) num2cell(entropyEMS)];

% Add to results
if ~exist('emptymean','var')
    result.EBE = cat(1,result.EBE,entropyE);
    result.ETE = cat(1,result.ETE,entropyES);
end
if ~exist('emptymedian','var')
    result.EBEM = cat(1,result.EBEM,entropyEM);
    result.ETEM = cat(1,result.ETEM,entropyEMS);
end

% ----------------------------
% Calculate BandPower Mappings
% ----------------------------
if exist('Mappings','var') & ~isempty(Mappings)
    for i = 1:size(Mappings.mappings,2);
        if size(Mappings.mappings{1,i},2) > 1;
            SpectMappingsMean(i,:)= mean(SpectAllMean(Mappings.mappings{1,i},:).^0.5).^2; %#ok<AGROW>
            SpectMappingsMedian(i,:)= median(SpectAllMedian(Mappings.mappings{1,i},:)); %#ok<AGROW>
        else
            SpectMappingsMean(i,:)= SpectAllMean(Mappings.mappings{1,i},:); %#ok<AGROW>
            SpectMappingsMedian(i,:)= SpectAllMedian(Mappings.mappings{1,i},:); %#ok<AGROW>
        end
    end
    % caluclate entropy
    entropy = zeros(size(Mappings.mappings,2),size(cfg.CollectFFT.freqs,1));
    entropyM = zeros(size(Mappings.mappings,2),size(cfg.CollectFFT.freqs,1));
    for i = 1:size(cfg.CollectFFT.freqs,1);
        range = cfg.CollectFFT.freqs(i,1):.01:cfg.CollectFFT.freqs(i,2);
        for j = 1:size(SpectMappingsMean,1)
            int_spec = exp(interp1(log(SpectAllF),log(SpectMappingsMean(j,:)),log(range),'linear','extrap'));
            entropy(j,i) = lab_spectral_entropy(int_spec);
            int_spec = exp(interp1(log(SpectAllF),log(SpectMappingsMedian(j,:)),log(range),'linear','extrap'));
            entropyM(j,i) = lab_spectral_entropy(int_spec);
        end
    end
    clearvars i j range int_spec
    
    % calculate total power
    range = min(cfg.CollectFFT.freqs(:,1)):.01:max(cfg.CollectFFT.freqs(:,2));
    for j = 1:size(SpectMappingsMean,1)
        int_spec = exp(interp1(log(SpectAllF),log(SpectMappingsMean(j,:)),log(range),'linear','extrap'));
        entropyS(j,1) = lab_spectral_entropy(int_spec); %#ok<AGROW>
        int_spec = exp(interp1(log(SpectAllF),log(SpectMappingsMedian(j,:)),log(range),'linear','extrap'));
        entropyMS(j,1) = lab_spectral_entropy(int_spec);  %#ok<AGROW>
    end
    clearvars j range int_spec
    
    % reshape to output format
    entropy = reshape(entropy,1,size(entropy,1)*size(entropy,2));
    entropyM = reshape(entropyM,1,size(entropyM,1)*size(entropyM,2));
    
    % add to results
    if ~exist('emptymean','var')
        entropy=[cellstr(calc.patient) num2cell(entropy)];
        entropyS=[cellstr(calc.patient) num2cell(entropyS')];
        result.BE = cat(1,result.BE,entropy);
        result.TE = cat(1,result.TE,entropyS);
    end
    if ~exist('emptymedian','var')
        entropyM=[cellstr(calc.patient) num2cell(entropyM)];
        entropyMS=[cellstr(calc.patient) num2cell(entropyMS')];
        result.BEM = cat(1,result.BEM,entropyM);
        result.TEM = cat(1,result.TEM,entropyMS);
    end
end

end