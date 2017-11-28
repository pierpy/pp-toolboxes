% Calculate Generalized linear model
%
% [Result,cfg] = lab_glm(data,outcome,header,cfg)
%
% data = array (subjects x measures)
% outcome = array (1 x subjects)
% header = see lab_read_statistics
% cfg (optional)
%
% Written by F. Hatz 2014

function [Result,cfg] = lab_glm(data,outcome,header,cfg)

if ~exist('cfg','var') | ~isfield(cfg,'GLM') | ~isfield(cfg.GLM,'distribution')
    if exist('outcome','var') & ~isempty(outcome) & islogical(outcome)
        cfg.GLM.distribution = 'binomial';
        cfg.GLM.link = 'logit';
    else
        cfg.GLM.distribution = 'normal';
        cfg.GLM.link = 'identity';
    end
end

if ~exist('data','var')
    % read data
    [data,header,result,~,cfg]  = lab_read_statistics(cfg,1,0,0,1,1);
    if isempty(data)
        Result = [];
        return
    end
    
    % settings
    [cfg,skipprocessing] = lab_set_glm(cfg,header,result);
    if skipprocessing == 1
        Result = [];
        return
    else
        pause(0.2);
    end
    
    numdata = size(data,2);
    Select = cfg.GLM.select(cfg.GLM.select<=numdata);
    if ~isempty(Select);
        data = data(:,Select);
        header.vars = header.vars(Select);
    end
    Select = cfg.GLM.select(cfg.GLM.select>numdata);
    Select = setdiff(Select-numdata,cfg.GLM.outcome);
    if ~isempty(Select);
        data = [data result(:,Select)];
        header.vars = [header.vars header.result(Select)];
    end
    outcome = result(:,cfg.GLM.outcome);
    header.outcome = header.result(cfg.GLM.outcome);
    tmp = find(isnan(outcome),1);
    if ~isempty(tmp)
        tmp = (~isnan(outcome));
        outcome = outcome(tmp);
        data = data(tmp,:);
        header.subjects = header.subjects(tmp);
    end
    clearvars Prompt Formats tmp
end
if exist('header','var') & isfield(header,'vars')
    varnames = header.vars;
    varnames2 = [varnames header.outcome];
else
    varnames = cellstr(num2str((1:size(data,2))'))';
    varnames2 = [varnames {'outcome'}];
end
if isfield(cfg,'file') & isfield(cfg,'path') & ~isempty(cfg.file) & ~isempty(cfg.path)
    [~,~,~,fileout] = lab_filename(cfg.file);
    warning off %#ok<WNOFF>
    mkdir(fullfile(cfg.path,cfg.GLM.folder));
    warning off %#ok<WNOFF>
    fileout = fullfile(fullfile(cfg.path,cfg.GLM.folder),fileout);
    cd(fullfile(cfg.path,cfg.GLM.folder));
else
    fileout = [];
end

% do calculation
if strcmp(cfg.GLM.distribution,'binomial')
    tmp = unique(outcome);
    outcome2 = NaN(size(outcome));
    outcome2(outcome==tmp(1)) = 0;
    outcome2(outcome==tmp(2)) = 1;
    Idx = find(~isnan(outcome2));
    outcome2 = outcome2(Idx,1);
    data = data(Idx,:);
    if exist('header','var') & isfield(header,'subjects')
        header.subjects = header.subjects(Idx,1);
    end
    outcome = logical(outcome2);
    cfg.islogical = true;
    clearvars outcome2 Idx
else
    cfg.islogical = false;
end

if strcmp(cfg.GLM.method,'stepwise')
    disp('Calculate GLM with stepwise reduction')
    if isfield(cfg.GLM,'bootstrap') & cfg.GLM.bootstrap > 0
        disp('Disable bootstrapping (only for fitting)')
        cfg.GLM.bootstrap = 0;
    end
    Mstats = GeneralizedLinearModel.stepwise(data,outcome,cfg.GLM.model, ...
        'Distribution',cfg.GLM.distribution,'Link',cfg.GLM.link, ...
        'VarNames',varnames2);
    stats.beta = Mstats.Coefficients.Estimate;
    stats.se = Mstats.Coefficients.SE;
    stats.t = Mstats.Coefficients.tStat;
    stats.p = Mstats.Coefficients.pValue;
    stats.dev = Mstats.Deviance;
    stats.dfe = Mstats.DFE;
    stats.sfit = [];
    stats.rsquared = Mstats.Rsquared;
    stats.criterion = Mstats.ModelCriterion;
    if length(Mstats.CoefficientNames) > 1
        varnames = Mstats.CoefficientNames(2:end);
    else
        varnames = [];
    end
    if Mstats.NumPredictors > 0
        D = devianceTest(Mstats);
        stats.D = dataset2struct(D);
        data = single(Mstats.Variables);
        data = data(:,1:end-1);
    else
        data = [];
    end
    if islogical(outcome) & ~isempty(data)
        Predict = predict(Mstats,data);
        if ~isfield(cfg.GLM,'ROC') | isempty(cfg.GLM.ROC)
            cfg.GLM.ROC.method = '>=';
        end
        cfg.GLM.ROC.filename = fileout;
        [stats.ROC,cfg.GLM.ROC] = lab_create_ROC(Predict,outcome,cfg.GLM.ROC);
    end
    stats.ModelStepwise = Mstats;
else
    disp('Calculate GLM')
    Mstats = GeneralizedLinearModel.fit(data,outcome, ...
        'Distribution',cfg.GLM.distribution,'Link',cfg.GLM.link, ...
        'VarNames',varnames2);
    stats.beta = Mstats.Coefficients.Estimate;
    stats.se = Mstats.Coefficients.SE;
    stats.t = Mstats.Coefficients.tStat;
    stats.p = Mstats.Coefficients.pValue;
    stats.dev = Mstats.Deviance;
    stats.dfe = Mstats.DFE;
    stats.sfit = [];
    stats.rsquared = Mstats.Rsquared;
    stats.criterion = Mstats.ModelCriterion;
    if length(Mstats.CoefficientNames) > 1
        varnames = Mstats.CoefficientNames(2:end);
    else
        varnames = [];
    end
    if Mstats.NumPredictors > 0
        D = devianceTest(Mstats);
        stats.D = dataset2struct(D);
        data = single(Mstats.Variables);
        data = data(:,1:end-1);
    else
        data = [];
    end
    
    %  [~,dev,stats] = glmfit(data,outcome,cfg.GLM.distribution,'link',cfg.GLM.link);
    %  stats.dev = dev;
    
    if islogical(outcome) & ~isempty(data)
        if size(data,2) > 1
            Predict = glmval(stats.beta,data,cfg.GLM.link);
        else
            disp('  Only one variable, use variable for calculation of ROC');
            Predict = data;
        end
        if ~isfield(cfg.GLM,'ROC') | isempty(cfg.GLM.ROC)
            cfg.GLM.ROC.method = '>=';
        end
        cfg.GLM.ROC.filename = fileout;
        [stats.ROC,cfg.GLM.ROC] = lab_create_ROC(Predict,outcome,cfg.GLM.ROC);
    end
    
    if isfield(cfg.GLM,'bootstrap') & cfg.GLM.bootstrap > 0 & ~isempty(data)
        disp('Calculate GLM bootstrapping')
        warning off %#ok<WNOFF>
        if isfield(cfg.GLM,'ROC') & isfield(cfg.GLM.ROC,'filename')
            cfg.GLM.ROC.filename = '';
        end
        [Stat_ci,Stat_B] = bootci(cfg.GLM.bootstrap,@do_glmfit,[data outcome]);
        warning on %#ok<WNON>
        Stat_B = reshape(Stat_B,[size(Stat_B,1) size(Stat_ci,2) size(Stat_ci,3)]);
        Stat_ci = permute(Stat_ci,[2 3 1]);
        Stat_B = permute(Stat_B,[2 3 1]);
        stats.Bootstrap.ci = Stat_ci;
        stats.Bootstrap.B = Stat_B;
        clearvars Stat_ci Stat_B
    end
end

if isfield(stats,'Bootstrap') & isfield(stats.Bootstrap,'B') & ~isempty(data)
    tmp = (stats.Bootstrap.B(:,4,:));
    tmp = reshape(tmp,size(tmp,1),size(tmp,3));
    labels = [{'constant'} regexprep(varnames,'_',' ')];
    fig1 = figure('NumberTitle','off','Name','GLM-Bootstrap-pValue','MenuBar','none','Color',[1 1 1]);
    m1 = uimenu(fig1,'Label','File');
    uimenu(m1,'Label','Save Picture','Callback',@(~,~)lab_print_figure);
    uimenu(m1,'Label','Copy Picture','Callback','print -dmeta;');
    uimenu(m1,'Label','Close','Callback','close;');
    m2 = uimenu(fig1,'Label','Edit');
    uimenu(m2,'Label','0 grad','Callback',@(~,~)set_X_Angle(labels,0));
    uimenu(m2,'Label','15 grad','Callback',@(~,~)set_X_Angle(labels,15));
    uimenu(m2,'Label','30 grad','Callback',@(~,~)set_X_Angle(labels,30));
    uimenu(m2,'Label','45 grad','Callback',@(~,~)set_X_Angle(labels,45));
    uimenu(m2,'Label','60 grad','Callback',@(~,~)set_X_Angle(labels,60));
    uimenu(m2,'Label','90 grad','Callback',@(~,~)set_X_Angle(labels,90));
    boxplot(tmp','labels',labels,'Colors','k','Symbol','+k');
    set_X_Angle(labels,45);
    if ~isempty(fileout)
        lab_print_figure([fileout '_Bootstrap_pValue'],fig1);
    end
    
    tmp = (stats.Bootstrap.B(:,1,:));
    tmp = reshape(tmp,size(tmp,1),size(tmp,3));
    labels = [{'constant'} regexprep(varnames,'_',' ')];
    fig1 = figure('NumberTitle','off','Name','GLM-Bootstrap-Beta','MenuBar','none','Color',[1 1 1]);
    m1 = uimenu(fig1,'Label','File');
    uimenu(m1,'Label','Save Picture','Callback',@(~,~)lab_print_figure);
    uimenu(m1,'Label','Copy Picture','Callback','print -dmeta;');
    uimenu(m1,'Label','Close','Callback','close;');
    m2 = uimenu(fig1,'Label','Edit');
    uimenu(m2,'Label','0 grad','Callback',@(~,~)set_X_Angle(labels,0));
    uimenu(m2,'Label','15 grad','Callback',@(~,~)set_X_Angle(labels,15));
    uimenu(m2,'Label','30 grad','Callback',@(~,~)set_X_Angle(labels,30));
    uimenu(m2,'Label','45 grad','Callback',@(~,~)set_X_Angle(labels,45));
    uimenu(m2,'Label','60 grad','Callback',@(~,~)set_X_Angle(labels,60));
    uimenu(m2,'Label','90 grad','Callback',@(~,~)set_X_Angle(labels,90));
    boxplot(tmp','labels',labels,'Colors','k','Symbol','+k');
    set_X_Angle(labels,45);
    if ~isempty(fileout)
        lab_print_figure([fileout '_Bootstrap_Beta'],fig1);
    end
    
    tmp = (stats.Bootstrap.B(1,5:end,:));
    tmp = reshape(tmp,size(tmp,2),size(tmp,3));
    fig1 = figure('NumberTitle','off','Name','GLM-Bootstrap','MenuBar','none','Color',[1 1 1]);
    m1 = uimenu(fig1,'Label','File');
    uimenu(m1,'Label','Save Picture','Callback',@(~,~)lab_print_figure);
    uimenu(m1,'Label','Copy Picture','Callback','print -dmeta;');
    uimenu(m1,'Label','Close','Callback','close;');
    m2 = uimenu(fig1,'Label','Edit');
    labels = {'AUC','Sensitivity','Spezificity','Threshold'};
    uimenu(m2,'Label','0 grad','Callback',@(~,~)set_X_Angle(labels,0));
    uimenu(m2,'Label','15 grad','Callback',@(~,~)set_X_Angle(labels,15));
    uimenu(m2,'Label','30 grad','Callback',@(~,~)set_X_Angle(labels,30));
    uimenu(m2,'Label','45 grad','Callback',@(~,~)set_X_Angle(labels,45));
    uimenu(m2,'Label','60 grad','Callback',@(~,~)set_X_Angle(labels,60));
    uimenu(m2,'Label','90 grad','Callback',@(~,~)set_X_Angle(labels,90));
    boxplot(tmp','labels',labels,'Colors','k','Symbol','+k');
    set_X_Angle(labels,45);
    if ~isempty(fileout)
        lab_print_figure([fileout '_Bootstrap'],fig1);
    end
    clearvars labels tmp tmp2
end

if ~isempty(fileout)
    % write result as xls-file
    resulttmp(:,1) = stats.beta;
    if isfield(stats,'se')
        resulttmp(:,2) = stats.se;
        resulttmp(:,3) = stats.t;
        resulttmp(:,4) = stats.p;
        xlsout = {'','B','SE','T','p-value'};
    else
        resulttmp(:,2) = stats.t;
        resulttmp(:,3) = stats.p;
        xlsout = {'','B','T','p-value'};
    end
    xlsout{2,1} = 'constant';
    if ~isempty(varnames)
        xlsout(3:2+length(varnames),1) = varnames';
    end
    xlsout(2:end,2:end) = num2cell(resulttmp);
    xlsout(end+1,1) = {''};
    xlsout(end+1,1:2) = {'degree of freedom',stats.dfe};
    if ~isempty(stats.sfit)
        xlsout(end+1,1:2) = {'sfit',stats.sfit};
    end
    if isfield(stats,'criterion')
        xlsout(end+1,1) = {''};
        xlsout(end+1,1) = {'-Criterion-'};
        xlsout(end+1,1:4) = {'AIC',stats.criterion.AIC,'corrected for sample size',stats.criterion.AICc};
        xlsout(end+1,1:2) = {'CAIC',stats.criterion.CAIC};
        xlsout(end+1,1:2) = {'BIC',stats.criterion.BIC};
    end
    if isfield(stats,'rsquared')
        xlsout(end+1,1) = {''};
        xlsout(end+1,1) = {'-R squared-'};
        xlsout(end+1,1:2) = {'R2 (ordinary)',stats.rsquared.Ordinary};
        xlsout(end+1,1:2) = {'R2 (adjusted)',stats.rsquared.Adjusted};
    end
    if isfield(stats,'D')
        if isfield(stats.D,'chi2Stat')
            xlsout(end+1,1) = {''};
            xlsout(end+1,1) = {'-Statistics-'};
            xlsout(end+1,1:2) = {'chi2stat',stats.D(end).chi2Stat};
            xlsout(end+1,1:2) = {'p-Value',stats.D(end).pValue};
        elseif isfield(stats.D,'FStat')
            xlsout(end+1,1) = {''};
            xlsout(end+1,1) = {'-Statistics-'};
            xlsout(end+1,1:2) = {'FStat',stats.D(end).FStat};
            xlsout(end+1,1:2) = {'p-Value',stats.D(end).pValue};
        end
    end
    if isfield(stats,'ROC') & isfield(stats.ROC,'AUC')
        xlsout(end+1,1) = {''};
        xlsout(end+1,1) = {'-ROC-'};
        xlsout(end+1,1:2) = {'AUC',stats.ROC.AUC};
        xlsout(end+1,1:2) = {'Sensitivity',stats.ROC.MaxYouden.Sensitivity};
        xlsout(end+1,1:2) = {'Specificity',stats.ROC.MaxYouden.Specificity};
        xlsout(end+1,1:2) = {'PPV',stats.ROC.MaxYouden.PPV};
        xlsout(end+1,1:2) = {'NPV',stats.ROC.MaxYouden.NPV};
        xlsout(end+1,1:2) = {'Threshold',stats.ROC.MaxYouden.Threshold};
    end
    filenameout = [fileout '_GLM.xls'];
    lab_write_xls(filenameout,xlsout);
    
    % write Matlab container-file
    filenameout = [fileout '_GLM.mat'];
    save(filenameout,'stats');
end

Result = stats;

    function Result = do_glmfit(data)
        outcomes = data(:,end);
        data = data(:,1:end-1);
        if length(unique(outcomes)) == 2 & min(outcomes) == 0 & max(outcomes) == 1
            outcomes = logical(outcomes);
        end
        Result = zeros(size(data,2)+1,8);
        % [~,~,statstmp] = glmfit(data,outcomes,cfg.GLM.distribution,'link',cfg.GLM.link);
        MstatsB = GeneralizedLinearModel.fit(data,outcomes, ...
            'Distribution',cfg.GLM.distribution,'Link',cfg.GLM.link, ...
            'VarNames',varnames2);
        statstmp.beta = MstatsB.Coefficients.Estimate;
        statstmp.se = MstatsB.Coefficients.SE;
        statstmp.t = MstatsB.Coefficients.tStat;
        statstmp.p = MstatsB.Coefficients.pValue;
        statstmp.dev = MstatsB.Deviance;
        statstmp.dfe = MstatsB.DFE;
        statstmp.sfit = [];
        statstmp.rsquared = MstatsB.Rsquared;
        statstmp.criterion = MstatsB.ModelCriterion;
        
        if islogical(outcomes) & ~isempty(data)
            if size(data,2) > 1
                Pred = glmval(statstmp.beta,data,cfg.GLM.link);
            else
                Pred = data;
            end
            if ~isfield(cfg.GLM,'ROC') | isempty(cfg.GLM.ROC)
                cfg.GLM.ROC.method = '>=';
            end
            [statstmp.ROC,cfg.GLM.ROC] = lab_create_ROC(Pred,outcomes,cfg.GLM.ROC);
        end
        
        Result(:,1) = statstmp.beta;
        if isfield(statstmp,'se')
            Result(:,2) = statstmp.se;
        end
        Result(:,3) = statstmp.t;
        Result(:,4) = statstmp.p;
        if isfield(statstmp,'ROC') & isfield(statstmp.ROC,'AUC')
            Result(1,5) = statstmp.ROC.AUC;
            Result(1,6) = statstmp.ROC.MaxYouden.Sensitivity;
            Result(1,7) = statstmp.ROC.MaxYouden.Specificity;
            Result(1,8) = statstmp.ROC.MaxYouden.Threshold;
        end
    end

end

function set_X_Angle(labels,angle)
   text_h = findobj(gca,'Type','text');
   if length(text_h) < length(labels)
       return
   end
   for cnt = 1:length(labels)
       if angle == 0
           set(text_h(length(text_h)-cnt+1),'Rotation',0, ...
               'String',labels{cnt},'HorizontalAlignment','center');
       elseif angle > 0
           set(text_h(length(text_h)-cnt+1),'Rotation',angle, ...
               'String',labels{cnt},'HorizontalAlignment','right');
       elseif angle < 0
           set(text_h(length(text_h)-cnt+1),'Rotation',angle, ...
               'String',labels{cnt},'HorizontalAlignment','left');
       end
   end
end