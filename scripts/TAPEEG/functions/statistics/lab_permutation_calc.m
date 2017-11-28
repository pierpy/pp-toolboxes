% Function calculates Permutation-Test
%
% [Rstat] = lab_permutation_calc(dataAll,resultAll,cfg,factorsAll)
%
%   dataAll   = Matrix (subjects x variables)
%   resultAll = vector (subjects x result)
%   cfg.permutations = number > 0
%   cfg.method = 'T-test' (unpaired t-test)
%                'Anova' (Anova)
%                'Spearman' (Spearman-rank-test)
%                'Kendall' (Kendall-Tau-rank-test)
%                'Pearson' (Pearson-correlation-test)
%                'Mann-Whitney-U' (Mann-Whitney-U-Test)
%                'Paired T-test' (paired t-test)
%   cfg.clustervars = number of variables in one cluster
%   factorsAll = factors to correct for by linear regression
%
% Written by F. Hatz 2011, Neurology University Hospital Basel

function [Rstat,skipprocessing] = lab_permutation_calc(dataAll,resultAll,cfg,factorsAll)

skipprocessing = 0;
Rstat = [];
if ~exist('cfg','var') | ~isfield(cfg,'permutations')
    cfg.permutations = 10000;
end

if size(dataAll,1) ~= size(resultAll,1)
    Rstat = [];
    skipprocessing = 1;
    disp 'Input data not valid'
    return
end

for Nresult = 1:size(resultAll,2)
    if length(cfg.permutations) >= Nresult
        permutations = cfg.permutations(Nresult);
    end
    if exist('factorsAll','var')
        factorNr = size(factorsAll,2);
        tmp = cat(2,dataAll,factorsAll,resultAll(:,Nresult));
    else
        factorNr = 0;
        tmp = cat(2,dataAll,resultAll(:,Nresult));
    end
    if ~strcmp(cfg.method{Nresult,1},'Paired T-test') & ~strcmp(cfg.method{Nresult,1},'Wilcoxon')
        tmp = sortrows(tmp,size(tmp,2));
    else
        [tmp1,tmp2,tmp3] = unique(tmp(:,end),'last');
        if length(tmp1) == 2 & (tmp2(2)-tmp2(1)) == 1
            tmp = cat(1,tmp(tmp3 == tmp1(1),:),tmp(tmp3 == tmp1(2),:));
            clearvars tmp1 tmp2 tmp3
        elseif strcmp(cfg.method{Nresult,1},'Paired T-test')
            disp('Permutation for paired t-test not possible, simple t-test is used')
            cfg.method{Nresult,1} = 'T-test';
            tmp = sortrows(tmp,size(tmp,2));
        else
            disp('Permutation for Wilcoxon not possible, Mann-Whitney-U is used')
            cfg.method{Nresult,1} = 'Mann-Whitney-U';
            tmp = sortrows(tmp,size(tmp,2));
        end
    end
    data = tmp(:,1:end-(factorNr+1));
    if factorNr  > 0
        factors = tmp(:,end-factorNr:end-1);
    end
    result = tmp(:,end);
    clearvars tmp
        
    tmp = intersect(find(~isnan(result)),find(min(~isnan(data),[],2)));
    data = data(tmp,:);
    result = result(tmp,:);
    if factorNr  > 0
        factors = factors(tmp,:);
    end
    [tmp,permctrlstep] = unique(result,'last');
    if ~isfield(cfg,'method')
        if length(permctrlstep) == 2
            cfg.method{Nresult,1} = 'T-test';
        elseif length(permctrlstep) > (size(data,1) / 3)
            cfg.method{Nresult,1} = 'Kendall';
        else
            cfg.method{Nresult,1} = 'Anova';
        end
    elseif strcmp(cfg.method,'T-test') & length(permctrlstep) > 2
        selection = listdlg('PromptString','Select results for T-test','SelectionMode','multiple','ListString',num2str(tmp(:)),'InitialValue',[1 2]);
        if length(selection) == 2
            selection = union(find(result==tmp(selection(1))),find(result==tmp(selection(2))));
            data = data(selection,:);
            if factorNr  > 0
                factors = factors(selection,:);
            end
            result = result(selection,:);
            [tmp,permctrlstep] = unique(result,'last');
        else
            disp('Error, less or more than 2 results selected, set method to ANOVA')
            cfg.method{Nresult,1} = 'Anova';
        end
    elseif strcmp(cfg.method,'Mann-Whitney-U') & length(permctrlstep) > 2
        selection = listdlg('PromptString','Select results for Mann-Whitney-U','SelectionMode','multiple','ListString',num2str(tmp(:)),'InitialValue',[1 2]);
        if length(selection) == 2
            selection = union(find(result==tmp(selection(1))),find(result==tmp(selection(2))));
            data = data(selection,:);
            if factorNr  > 0
                factors = factors(selection,:);
            end
            result = result(selection,:);
            [tmp,permctrlstep] = unique(result,'last');
        else
            disp('Error, less or more than 2 results selected, set method to KruskalWallis')
            cfg.method{Nresult,1} = 'KruskalWallis';
        end
    end
    clearvars tmp
    if permutations > 1
        disp (['Permutation with ' cfg.method{Nresult,1} ', ' num2str(permutations) ' permutations'])
    end
    
    % Do linear regression to correct for factors
    if factorNr > 0
        [dataR,stat.stats] = lab_linearregression(data,factors,result);
    else
        dataR = data;
        stat.stats = [];
    end
    
    rng('default');
    rng('shuffle');
    permutationstmp = floor(permutations/20);
    if strcmp(cfg.method{Nresult,1},'Paired T-test')
        if  size(result,1)/2 == floor(size(result,1)/2) & permctrlstep(1) == size(result,1)/2
            fprintf(' Calculate permutation for paired t-test (only paired subjects are flipped)')
            % Insert baseline as first permutation
            permcalc(:,1) = 1:size(result,1);
            permdouble = 0;
            for i = 2:permutations
                k = 1;
                while k == 1
                    selection = rand(size(result,1)/2,1);
                    selection = find(selection > 0.5);
                    tmp = permcalc(:,1);
                    tmp(selection,1) = selection + size(result,1)/2;
                    tmp(selection + size(result,1)/2,1) = selection;
                    k = 0;
                    for j = 1:(i-1)
                        if isequal(tmp(:,1),permcalc(:,j))
                            k = 1;
                            permdouble = permdouble + 1;
                            if permdouble > permutations*2
                                skipprocessing = 1;
                                disp('Error: Permutation not possible (you need more subjects)')
                                return
                            end
                        end
                    end
                end
                if mod(i,permutationstmp) == 0
                    fprintf('.')
                end
                permcalc(:,i) = tmp(:,1); %#ok<AGROW>
            end
            clearvars 'tmp' 'i' 'j' 'calc' 'k' 'selection';
            disp(':')
            disp ('Calculate permutation matrix finished')
        else
            disp('Permutation for paired t-test not possible, simple t-test is used')
            cfg.method{Nresult,1} = 'T-test';
        end
    end
    
    if permutations > 1 & (~exist('permcalc','var') | size(permcalc,1) ~= size(result,1))
        fprintf('Calculate permutation matrix')
        % prepare matrix for permutation
        calc = result;
        calc(:,2) = 1:size(calc,1);
        permcalc = zeros(size(calc,1),permutations);
        % Insert baseline as first permutation
        permcalc(:,1) = calc(:,2);
        % Build permutation matrix and control/exclude repeated permutations
        permdouble = 0;
        for i = 2:permutations
            k = 1;
            while k == 1
                calc(:,3) = rand(size(calc,1),1);
                tmp = sortrows(calc,3);
                tmp(:,1) = calc(:,1);
                tmp = sortrows(tmp,[1 2]);
                k = 0;
                for j = 1:(i-1)
                    if isequal(tmp(:,2),permcalc(:,j))
                        k = 1;
                        permdouble = permdouble + 1;
                        if permdouble > permutations
                            disp(':')
                            disp('Error: Permutation not possible (you need more subjects)')
                            skipprocessing = 1;
                            return
                        end
                    end
                end
            end
            if mod(i,permutationstmp) == 0
                fprintf('.')
            end
            permcalc(:,i) = tmp(:,2);
        end
        clearvars 'tmp' 'i' 'j' 'calc' 'k';
        disp(':')
        disp ('Calculate permutation matrix finished')
    elseif permutations == 1
        permcalc = (1:size(result,1))';
    end
    
    %---------------------
    % Calculate statistics
    %---------------------
    perm.p = zeros(permutations,size(dataR,2));
    
    % Calculate t-test
    if strcmp(cfg.method{Nresult,1},'T-test')
        fprintf('  T-test')
        perm.T = zeros(permutations,size(dataR,2));
        for i = 1:permutations
            tmp = dataR(permcalc(:,i),:);
            [~,p,~,stats] = ttest2(tmp(1:permctrlstep(1),:),tmp((permctrlstep(1)+1):end,:));
            perm.T(i,:) = stats.tstat;
            perm.p(i,:) = p;
            perm.df = stats.df;
            if mod(i,permutationstmp) == 0
                fprintf('.')
            end
        end
        disp(':')
        clearvars 'h' 'p' 'ci' 'stats' 'i'
        stat.T = perm.T(1,:);
        stat.p = perm.p(1,:);
        stat.df = perm.df(1,:);
        stat = store_mean_std(dataR,permctrlstep,stat,result);
    end
    
    % Calculate paired t-test
    if strcmp(cfg.method{Nresult,1},'Paired T-test')
        fprintf('  Paired t-test')
        perm.T = zeros(permutations,size(dataR,2));
        for i = 1:permutations
            tmp = dataR(permcalc(:,i),:);
            [~,p,~,stats] = ttest(tmp(1:permctrlstep(1),:),tmp((permctrlstep(1)+1):end,:));
            perm.T(i,:) = stats.tstat;
            perm.p(i,:) = p;
            perm.df = stats.df;
            if mod(i,permutationstmp) == 0
                fprintf('.')
            end
        end
        disp(':')
        clearvars 'h' 'p' 'ci' 'stats' 'i'
        stat.T = perm.T(1,:);
        stat.p = perm.p(1,:);
        stat.df = perm.df(1,:);
        stat = store_mean_std(dataR,permctrlstep,stat,result);
    end
    
    % Calculate Mann-Whitney-U
    if strcmp(cfg.method{Nresult,1},'Mann-Whitney-U')
        fprintf('   Mann-Whitney-U')
        perm.Z = zeros(permutations,size(dataR,2));
        for i = 1:permutations
            tmp = dataR(permcalc(:,i),:);
            for j = 1:size(tmp,2)
                [p,~,stats] = ranksum(tmp(1:permctrlstep(1),j),tmp((permctrlstep(1)+1):end,j));
                perm.Z(i,j) = stats.zval;
                perm.p(i,j) = p;
            end
            if mod(i,permutationstmp) == 0
                fprintf('.')
            end
        end
        disp(':')
        clearvars 'h' 'p' 'ci' 'stats' 'i'
        stat.Z = perm.Z(1,:);
        stat.p = perm.p(1,:);
        stat = store_mean_std(dataR,permctrlstep,stat,result);
    end
    
    % Calculate Wilcoxon
    if strcmp(cfg.method{Nresult,1},'Wilcoxon')
        fprintf('   Wilcoxon')
        perm.Z = zeros(permutations,size(dataR,2));
        for i = 1:permutations
            tmp = dataR(permcalc(:,i),:);
            for j = 1:size(tmp,2)
                [p,~,stats] = signrank(tmp(1:permctrlstep(1),j),tmp((permctrlstep(1)+1):end,j));
                if isfield(stats,'zval')
                    perm.Z(i,j) = stats.zval;
                else
                    perm.Z(i,j) = 0;
                end
                perm.p(i,j) = p;
            end
            if mod(i,permutationstmp) == 0
                fprintf('.')
            end
        end
        disp(':')
        clearvars 'h' 'p' 'ci' 'stats' 'i'
        stat.Z = perm.Z(1,:);
        stat.p = perm.p(1,:);
        stat = store_mean_std(dataR,permctrlstep,stat,result);
    end
    
    % Calculate anova
    if strcmp(cfg.method{Nresult,1},'Anova')
        fprintf('  Anova')
        perm.F = zeros(permutations,size(dataR,2));
        for i = 1:permutations
            tmp = dataR(permcalc(:,i),:);
            for j = 1:size(tmp,2)
                [p,table] = anova1(tmp(:,j),result,'off');
                perm.F(i,j) = table{2,5};
                perm.p(i,j) = p;
            end
            if mod(i,permutationstmp) == 0
                fprintf('.')
            end
        end
        disp(':')
        clearvars 'table' 'p' 'i'
        stat.F = perm.F(1,:);
        stat.p = perm.p(1,:);
        stat = store_mean_std(dataR,permctrlstep,stat,result);
    end
    
    % Calculate RM-Anova
    if strcmp(cfg.method{Nresult,1},'RM-Anova')
        fprintf('  RM-Anova')
        if ~isfield(cfg,'clustervars') | isempty(cfg.clustervars) | cfg.clustervars > size(dataR,2)
            cfg.clustervars = size(dataR,2);
        end
        if ~isfield(cfg,'preservemean') | isempty(cfg.preservemean)
            cfg.preservemean = true;
        end
        NVar = floor(size(dataR,2) / cfg.clustervars);
        perm.F = zeros(permutations,NVar);
        Ridx = unique(result);
        if length(Ridx) > 1
            cfg.AnovaVars = {'Time','Group','Ineratcion','Subjects (matching)'};
        else
            cfg.AnovaVars = {'Time','Subjects (matching)'};
        end
        for i = 1:permutations
            datatmp = cell(1,length(Ridx));
            for j = 1:length(Ridx)
                for m = 1:NVar
                    if cfg.preservemean == true
                        % preserve mean value per subject/trial while doing permutation
                        tmp = dataR(:,(m-1)*cfg.clustervars+1:m*cfg.clustervars);
                        tmp2 = mean(tmp,2);
                        tmp = tmp - repmat(tmp2,1,size(tmp,2));
                        tmp = tmp(permcalc(:,i),:);
                        tmp = tmp + repmat(tmp2,1,size(tmp,2));
                        datatmp{m,j} = tmp(result == Ridx(j),:);
                    else
                        datatmp{m,j} = tmp(result == Ridx(j),(m-1)*cfg.clustervars+1:m*cfg.clustervars);
                    end
                end
            end
            for j = 1:NVar
                [~,table] = anova_rm(datatmp(j,:),'off');
                perm.F(i,j,1) = table{2,5};
                perm.F(i,j,1) = table{2,6};
                perm.F(i,j,2) = table{3,5};
                perm.F(i,j,2) = table{3,6};
                if size(table,1) > 5
                    perm.F(i,j,3) = table{4,5};
                    perm.F(i,j,3) = table{4,6};
                    perm.F(i,j,4) = table{5,5};
                    perm.F(i,j,4) = table{5,6};
                end
            end
            if mod(i,permutationstmp) == 0
                fprintf('.')
            end
        end
        disp(':')
        clearvars 'table' 'p' 'i'
        stat.F = perm.F(1,:);
        stat.p = perm.p(1,:);
        stat = store_mean_std(dataR,permctrlstep,stat,result);
    end
    
    % Calculate KruskalWallis
    if strcmp(cfg.method{Nresult,1},'KruskalWallis')
        fprintf('  KruskalWallis')
        perm.F = zeros(permutations,size(dataR,2));
        for i = 1:permutations
            tmp = dataR(permcalc(:,i),:);
            for j = 1:size(tmp,2)
                [p,table] = kruskalwallis(tmp(:,j),result,'off');
                perm.F(i,j) = table{2,5};
                perm.p(i,j) = p;
            end
            if mod(i,permutationstmp) == 0
                fprintf('.')
            end
        end
        disp(':')
        clearvars 'table' 'p' 'i'
        stat.F = perm.F(1,:);
        stat.p = perm.p(1,:);
        stat = store_mean_std(dataR,permctrlstep,stat,result);
    end
    
    % Calculate corr
    if strcmp(cfg.method{Nresult,1},'Kendall') || strcmp(cfg.method{Nresult,1},'Spearman') || strcmp(cfg.method{Nresult,1},'Pearson')
        fprintf(['  ' cfg.method{Nresult,1}])
        perm.R = zeros(permutations,size(dataR,2));
        for i = 1:permutations
            tmp = dataR(permcalc(:,i),:);
            [r,p] = corr(tmp,result,'type',cfg.method{Nresult,1});
            perm.R(i,:) = r;
            perm.p(i,:) = p;
            if mod(i,permutationstmp) == 0
                fprintf('.')
            end
        end
        disp(':')
        clearvars 'r' 'p' 'i'
        stat.R = perm.R(1,:);
        stat.p = perm.p(1,:);
    end
    
    % Calculate corrected pValue
    if permutations > 1
        if ~strcmp(cfg.method{Nresult,1},'RM-Anova')
            minp = min(perm.p,[],2);
            if minp == 0;
                error 'Permutation was not working (method?)'
            end
            minp = repmat(minp,1,size(perm.p,2));
            tmp = (minp < repmat(perm.p(1,:),permutations,1));
            tmp = sum(tmp,1);
            stat.mxp = tmp / size(perm.p,1);
            clearvars 'tmp' 'minp'
            if isfield(cfg,'clustervars2') & cfg.clustervars2 > 1
                stat.mxpV = zeros(1,size(perm.p,2));
                step = cfg.clustervars * cfg.clustervars2;
                for i = 1:(size(perm.p,2) / step)
                    tmp.start = (i-1)*step + 1;
                    tmp.end = i * step;
                    minp = min(perm.p(:,tmp.start:tmp.end),[],2);
                    minp = repmat(minp,1,step);
                    tmp.result=(minp < repmat(perm.p(1,tmp.start:tmp.end),permutations,1));
                    tmp.result = sum(tmp.result);
                    stat.mxpV(1,tmp.start:tmp.end) = tmp.result / size(perm.p,1);
                    clearvars tmp minp
                end
            end
            if isfield(cfg,'clustervars') & cfg.clustervars > 1
                stat.mxpC = zeros(1,size(perm.p,2));
                for i = 1:(size(perm.p,2) / cfg.clustervars)
                    tmp.start = (i-1)*cfg.clustervars + 1;
                    tmp.end = i * cfg.clustervars;
                    minp = min(perm.p(:,tmp.start:tmp.end),[],2);
                    minp = repmat(minp,1,cfg.clustervars);
                    tmp.result=(minp < repmat(perm.p(1,tmp.start:tmp.end),permutations,1));
                    tmp.result = sum(tmp.result);
                    stat.mxpC(1,tmp.start:tmp.end) = tmp.result / size(perm.p,1);
                    clearvars tmp minp
                end
            end
            stat.mxpS = zeros(1,size(perm.p,2));
            for i = 1:size(perm.p,2)
                result=(perm.p(:,i) < repmat(perm.p(1,i),permutations,1));
                result = sum(result);
                stat.mxpS(1,i) = result / size(perm.p,1);
                clearvars result
            end
        else
            stat.mxp = zeros(size(perm.p,3),size(perm.p,2));
            if isfield(cfg,'clustervars2') & cfg.clustervars2 > 1
                stat.mxpV = zeros(size(perm.p,3),size(perm.p,2));
            end
            stat.mxpS = zeros(size(perm.p,3),size(perm.p,2));
            for j = 1:size(perm.p,3)
                minp = min(perm.p(:,:,j),[],2);
                if minp == 0;
                    error 'Permutation was not working (method?)'
                end
                minp = repmat(minp,1,size(perm.p,2));
                tmp = (minp < repmat(perm.p(1,:,j),permutations,1));
                tmp = sum(tmp,1);
                stat.mxp(j,:) = tmp / size(perm.p,1);
                clearvars 'tmp' 'minp'
                if isfield(cfg,'clustervars2') & cfg.clustervars2 > 1
                    stat.mxpV = zeros(1,size(perm.p,2));
                    step = cfg.clustervars;
                    for i = 1:(size(perm.p,2) / step)
                        tmp.start = (i-1)*step + 1;
                        tmp.end = i * step;
                        minp = min(perm.p(:,tmp.start:tmp.end,j),[],2);
                        minp = repmat(minp,1,step);
                        tmp.result=(minp < repmat(perm.p(1,tmp.start:tmp.end,j),permutations,1));
                        tmp.result = sum(tmp.result);
                        stat.mxpV(j,tmp.start:tmp.end) = tmp.result / size(perm.p,1);
                        clearvars tmp minp
                    end
                end
                stat.mxpS = zeros(1,size(perm.p,2));
                for i = 1:size(perm.p,2)
                    result = (perm.p(:,i,j) < repmat(perm.p(1,i,j),permutations,1));
                    result = sum(result);
                    stat.mxpS(j,i) = result / size(perm.p,1);
                    clearvars result
                end
            end
            stat.mxp = stat.mxp';
            stat.mxp = stat.mxp(:)';
            if isfield(cfg,'clustervars2') & cfg.clustervars2 > 1
                stat.mxpV = stat.mxpV';
                stat.mxpV = stat.mxpV(:)';
            end
            stat.mxpS = stat.mxpS';
            stat.mxpS = stat.mxpS(:)';
        end
    end
    
    % Store results in Rstat
    if size(resultAll,2) == 1
        Rstat = stat;
    else
        Rstat{Nresult} = stat; %#ok<AGROW>
    end
    clearvars stat permcalc
end

    function stat = store_mean_std(dataR,permctrlstep,stat,result)
        for ndat = 1:size(dataR,2)
            tmpU = unique(dataR(~isnan(dataR(:,ndat)),ndat));
            tmpR = dataR(1:permctrlstep(1),ndat);
            tmpR = tmpR(~isnan(tmpR));
            if length(tmpU) == 2
                stat.logical{ndat} = tmpU(:)';
                stat.mean{1}(ndat) = sum(tmpR==stat.logical{ndat}(1));
                stat.std{1}(ndat) = sum(tmpR==stat.logical{ndat}(2));
            else
                stat.logical{ndat} = [];
                stat.mean{1}(ndat) = mean(tmpR);
                stat.std{1}(ndat) = std(tmpR);
                stat.median{1}(ndat) = median(tmpR);
                tmpR = sort(tmpR(:));
                stat.percent25{1}(ndat) = tmpR(ceil(0.25*length(tmpR)));
                stat.percent75{1}(ndat) = tmpR(ceil(0.75*length(tmpR)));
                stat.min{1}(ndat) = min(tmpR);
                stat.max{1}(ndat) = max(tmpR);
            end
        end
        stat.group{1} = num2str(result(1,1));
        for ndat = 1:size(dataR,2)
            for ngroup = 1:length(permctrlstep)-1
                tmpR = dataR(permctrlstep(ngroup)+1:permctrlstep(ngroup+1),ndat);
                tmpR = tmpR(~isnan(tmpR));
                if ~isempty(stat.logical{ndat})
                    stat.mean{ngroup+1}(ndat) = sum(tmpR==stat.logical{ndat}(1));
                    stat.std{ngroup+1}(ndat) = sum(tmpR==stat.logical{ndat}(2));
                else
                    stat.mean{ngroup+1}(1,ndat) = mean(tmpR);
                    stat.std{ngroup+1}(1,ndat) = std(tmpR);
                    stat.median{ngroup+1}(ndat) = median(tmpR);
                    tmpR = sort(tmpR(:));
                    stat.percent25{ngroup+1}(ndat) = tmpR(ceil(0.25*length(tmpR)));
                    stat.percent75{ngroup+1}(ndat) = tmpR(ceil(0.75*length(tmpR)));
                    stat.min{ngroup+1}(ndat) = min(tmpR);
                    stat.max{ngroup+1}(ndat) = max(tmpR);
                end
            end
        end
        for ngroup = 1:length(permctrlstep)-1
            stat.group{ngroup+1} = num2str(result(permctrlstep(ngroup)+1,1));
        end
    end

end

