% Helper function for lab_collect_spectras
%
% [cfg,result] = lab_collect_spectraldata_bp(cfg,result,calc,SpectAll,SpectAllMean,SpectAllMedian,SpectAllF,Mappings)
%
% written by F. Hatz 2012

function [cfg,result] = lab_collect_spectraldata_bp(cfg,result,calc,SpectAll,SpectAllMean,SpectAllMedian,SpectAllF,Mappings)

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
bandpowerE = zeros(size(SpectAll,3),size(cfg.CollectFFT.freqs,1));
bandpowerEM = zeros(size(SpectAll,3),size(cfg.CollectFFT.freqs,1));
for i = 1:(size(cfg.CollectFFT.freqs,1));
    range = cfg.CollectFFT.freqs(i,1):.01:cfg.CollectFFT.freqs(i,2);
    for j = 1:size(SpectAllMean,1)
        int_spec = exp(interp1(log(SpectAllF),log(SpectAllMean(j,:)),log(range),'linear','extrap'));
        bandpowerE(j,i) = trapz(range,int_spec);
        int_spec = exp(interp1(log(SpectAllF),log(SpectAllMedian(j,:)),log(range),'linear','extrap'));
        bandpowerEM(j,i) = trapz(range,int_spec);
    end
end
clearvars j range int_spec

% calculate total power
range = min(cfg.CollectFFT.freqs(:,1)):.01:max(cfg.CollectFFT.freqs(:,2));
for j = 1:size(SpectAllMean,1)
    int_spec = exp(interp1(log(SpectAllF),log(SpectAllMean(j,:)),log(range),'linear','extrap'));
    bandpowerES(j,1) = trapz(range,int_spec); %#ok<AGROW>
    int_spec = exp(interp1(log(SpectAllF),log(SpectAllMedian(j,:)),log(range),'linear','extrap'));
    bandpowerEMS(j,1) = trapz(range,int_spec); %#ok<AGROW>
end
clearvars j range int_spec

% calculate realitve power
bandpowerER = bandpowerE ./ repmat(bandpowerES,1,size(bandpowerE,2));
bandpowerEMR = bandpowerEM ./ repmat(bandpowerEMS,1,size(bandpowerEM,2));

% reshape data to output format
bandpowerE = reshape(bandpowerE,1,size(bandpowerE,1)*size(bandpowerE,2));
bandpowerEM = reshape(bandpowerEM,1,size(bandpowerEM,1)*size(bandpowerEM,2));
bandpowerEMR = reshape(bandpowerEMR,1,size(bandpowerEMR,1)*size(bandpowerEMR,2));
bandpowerER = reshape(bandpowerER,1,size(bandpowerER,1)*size(bandpowerER,2));

% log and logit data
bandpowerELog=[cellstr(calc.patient) num2cell(log(bandpowerE))];
bandpowerEMLog=[cellstr(calc.patient) num2cell(log(bandpowerEM))];
bandpowerEMRLog=[cellstr(calc.patient) num2cell(lab_logit(bandpowerEMR))];
bandpowerERLog=[cellstr(calc.patient) num2cell(lab_logit(bandpowerER))];
bandpowerE=[cellstr(calc.patient) num2cell(bandpowerE)];
bandpowerEM=[cellstr(calc.patient) num2cell(bandpowerEM)];
bandpowerEMR=[cellstr(calc.patient) num2cell(bandpowerEMR)];
bandpowerER=[cellstr(calc.patient) num2cell(bandpowerER)];

% Add to results
if ~exist('emptymean','var')
    result.EBP = cat(1,result.EBP,bandpowerE);
    result.EBPR = cat(1,result.EBPR,bandpowerER);
    result.EBPlog = cat(1,result.EBPlog,bandpowerELog);
    result.EBPRlog = cat(1,result.EBPRlog,bandpowerERLog);
end
if ~exist('emptymedian','var')
    result.EBPM = cat(1,result.EBPM,bandpowerEM);
    result.EBPMR = cat(1,result.EBPMR,bandpowerEMR);
    result.EBPMlog = cat(1,result.EBPMlog,bandpowerEMLog);
    result.EBPMRlog = cat(1,result.EBPMRlog,bandpowerEMRLog);
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
    % caluclate bandpower
    bandpower = zeros(size(Mappings.mappings,2),size(cfg.CollectFFT.freqs,1));
    bandpowerM = zeros(size(Mappings.mappings,2),size(cfg.CollectFFT.freqs,1));
    for i = 1:size(cfg.CollectFFT.freqs,1);
        range = cfg.CollectFFT.freqs(i,1):.01:cfg.CollectFFT.freqs(i,2);
        for j = 1:size(SpectMappingsMean,1)
            int_spec = exp(interp1(log(SpectAllF),log(SpectMappingsMean(j,:)),log(range),'linear','extrap'));
            bandpower(j,i) = trapz(range,int_spec);
            int_spec = exp(interp1(log(SpectAllF),log(SpectMappingsMedian(j,:)),log(range),'linear','extrap'));
            bandpowerM(j,i) = trapz(range,int_spec);
        end
    end
    clearvars i j range int_spec
    
    % calculate total power
    range = min(cfg.CollectFFT.freqs(:,1)):.01:max(cfg.CollectFFT.freqs(:,2));
    for j = 1:size(SpectMappingsMean,1)
        int_spec = exp(interp1(log(SpectAllF),log(SpectMappingsMean(j,:)),log(range),'linear','extrap'));
        bandpowerS(j,1) = trapz(range,int_spec); %#ok<AGROW>
        int_spec = exp(interp1(log(SpectAllF),log(SpectMappingsMedian(j,:)),log(range),'linear','extrap'));
        bandpowerMS(j,1) = trapz(range,int_spec); %#ok<AGROW>
    end
    clearvars j range int_spec
    
    % calculate realitve power
    bandpowerR = bandpower ./ repmat(bandpowerS,1,size(bandpower,2));
    bandpowerMR = bandpowerM ./ repmat(bandpowerMS,1,size(bandpowerM,2));
    
    % reshape to output format
    bandpower = reshape(bandpower,1,size(bandpower,1)*size(bandpower,2));
    bandpowerM = reshape(bandpowerM,1,size(bandpowerM,1)*size(bandpowerM,2));
    bandpowerMR = reshape(bandpowerMR,1,size(bandpowerMR,1)*size(bandpowerMR,2));
    bandpowerR = reshape(bandpowerR,1,size(bandpowerR,1)*size(bandpowerR,2));
    
    % log and logit data
    bandpowerLog=[cellstr(calc.patient) num2cell(log(bandpower))];
    bandpower=[cellstr(calc.patient) num2cell(bandpower)];
    bandpowerMLog=[cellstr(calc.patient) num2cell(log(bandpowerM))];
    bandpowerM=[cellstr(calc.patient) num2cell(bandpowerM)];
    bandpowerMRLog=[cellstr(calc.patient) num2cell(lab_logit(bandpowerMR))];
    bandpowerMR=[cellstr(calc.patient) num2cell(bandpowerMR)];
    bandpowerRLog=[cellstr(calc.patient) num2cell(lab_logit(bandpowerR))];
    bandpowerR=[cellstr(calc.patient) num2cell(bandpowerR)];
    
    % add to results
    if ~exist('emptymean','var')
        result.BP = cat(1,result.BP,bandpower);
        result.BPR = cat(1,result.BPR,bandpowerR);
        result.BPlog = cat(1,result.BPlog,bandpowerLog);
        result.BPRlog = cat(1,result.BPRlog,bandpowerRLog);
    end
    if ~exist('emptymedian','var')
        result.BPM = cat(1,result.BPM,bandpowerM);
        result.BPMR = cat(1,result.BPMR,bandpowerMR);
        result.BPMlog = cat(1,result.BPMlog,bandpowerMLog);
        result.BPMRlog = cat(1,result.BPMRlog,bandpowerMRLog);
    end
end

end