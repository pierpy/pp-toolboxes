function [Result,settings] = lab_create_ROC(data,outcome,settings)
    
    if ~exist('data','var')
        % read data
        [data,~,outcome,~,cfg]  = lab_read_statistics([],1,0,0,1,1);
        if isempty(data)
            Result = [];
            settings = [];
            return
        end
        if isfield(cfg,'file') & isfield(cfg,'path') & ~isempty(cfg.file) & ~isempty(cfg.path)
            Filename = fullfile(cfg.path,cfg.file);
        end
    end
    
    if ~exist('settings','var')
        settings = [];
        [settings,skipprocessing] = lab_set_create_ROC(settings);
        if skipprocessing == 1
            Result = [];
            return
        else
            pause(0.2);
        end
        if exist('Filename','var')
            settings.filename = Filename;
        end
    end
    
    outcome = outcome(:);
    if size(data,2) == size(outcome,1)
        data = data';
    end
    if size(data,2) > 1
        disp('  Input has more than one variable, reduce to first variable')
        data = data(:,1);
    end
    
    if ~isfield(settings,'method')
        settings.method = '<=';
    end
    if ~isnumeric(outcome) & ~islogical(outcome)
        [tmp,~,tmp2] = unique(outcome);
        if isfield(settings,'target') & ~isempty(settings.target)
            tmp3 = find(strcmp(tmp,settings.target) == 1);
            if isempty(tmp3)
                settings.target = 1;
            elseif min(tmp3) > 2
                Idx = union(find(strcmp(outcome,tmp(tmp3(1))) == 1),find(strcmp(outcome,tmp(1)) == 1));
                outcome = outcome(Idx,:);
                data = data(Idx,:);
                clearvars Idx
                [~,~,tmp2] = unique(outcome);
                settings.target = 1;
            else
                settings.target = tmp3(1) - 1;
            end
        else
            settings.target = 1;
        end
        tmp2 = tmp2 - 1;
        Idx = find(tmp2 <= 1);
        data = data(Idx);
        outcome = logical(tmp2(Idx));
    elseif ~islogical(outcome)
        [tmp,~,tmp2] = unique(outcome);
        if isfield(settings,'target') & ~isempty(settings.target)
            tmp3 = find(tmp == settings.target);
            if isempty(tmp3)
                settings.target = 1;
            elseif min(tmp3) > 2
                Idx = union(find(outcome == tmp(tmp3(1))),find(outcome == tmp(1)));
                outcome = outcome(Idx,:);
                data = data(Idx,:);
                clearvars Idx
                [~,~,tmp2] = unique(outcome);
                settings.target = 1;
            else
                settings.target = tmp3(1) - 1;
            end
        else
            settings.target = 1;
        end
        tmp2 = tmp2 - 1;
        Idx = find(tmp2 <= 1);
        data = data(Idx);
        outcome = logical(tmp2(Idx));
    elseif isfield(settings,'target') & settings.target == false
        outcome = logical(abs(outcome-1));
    else
        settings.target = 1;
    end
    
    Ndata = length(outcome);
    Thresholds = sort(data,'descend');
    Resulttmp = zeros(Ndata,8);
    for i = 1:Ndata
        Predict = false(Ndata,1);
        if strcmp(settings.method,'<=')
            Predict(data <= Thresholds(i)) = true;
        elseif strcmp(settings.method,'<')
            Predict(data < Thresholds(i)) = true;
        elseif strcmp(settings.method,'>=')
            Predict(data >= Thresholds(i)) = true;
        else
            Predict(data > Thresholds(i)) = true;
        end
        Resulttmp(i,:) = ROCevaluate([outcome Predict]);
    end
    Result.Accuracy = Resulttmp(:,1)';
    Result.Sensitivity = Resulttmp(:,2)';
    Result.Specificity = Resulttmp(:,3)';
    Result.PPV = Resulttmp(:,4)';
    Result.NPV = Resulttmp(:,5)';
    Result.Precision = Resulttmp(:,6)';
    Result.F = Resulttmp(:,7)';
    Result.Gmean = Resulttmp(:,8)';
    if strcmp(settings.method,'<=') | strcmp(settings.method,'<')
        [~,~,~,Result.AUC] = perfcurve(outcome,-data,true);
    else
        [~,~,~,Result.AUC] = perfcurve(outcome,data,true);
    end
    
    Result.Youden = Result.Sensitivity + Result.Specificity;
    MaxYouden = find(Result.Youden == max(Result.Youden));
    Result.MaxYouden.Threshold = Thresholds(MaxYouden(1));
    Result.MaxYouden.Sensitivity = Result.Sensitivity(MaxYouden(1));
    Result.MaxYouden.Specificity = Result.Specificity(MaxYouden(1));
    Result.MaxYouden.PPV = Result.PPV(MaxYouden(1));
    Result.MaxYouden.NPV = Result.NPV(MaxYouden(1));
    if isfield(settings,'nbootstrap') & settings.nbootstrap > 0
        Predict = false(Ndata,1);
        if strcmp(settings.method,'<=')
            Predict(data <= Result.Threshold) = true;
        elseif strcmp(settings.method,'<')
            Predict(data < Result.Threshold) = true;
        elseif strcmp(settings.method,'>=')
            Predict(data >= Result.Threshold) = true;
        else
            Predict(data > Result.Threshold) = true;
        end
        BBci = bootci(settings.nbootstrap,@ROCevaluate,[outcome Predict]);
        Result.BOOTSTRAP.Accuracy = BBci(:,1)';
        Result.BOOTSTRAP.Sensitivity = BBci(:,2)';
        Result.BOOTSTRAP.Specificity = BBci(:,3)';
        Result.BOOTSTRAP.PPV = BBci(:,4)';
        Result.BOOTSTRAP.NPV = BBci(:,5)';
        Result.BOOTSTRAP.Precision = BBci(:,6)';
        Result.BOOTSTRAP.F = BBci(:,7)';
        Result.BOOTSTRAP.Gmean = BBci(:,8)';
        clearvars BBci
        if strcmp(settings.method,'<=') | strcmp(settings.method,'<')
            BBci = bootci(settings.nbootstrap,@AUCevaluate,[outcome -data]);
        else
            BBci = bootci(settings.nbootstrap,@AUCevaluate,[outcome data]);
        end
        Result.BOOTSTRAP.AUC = BBci';
    end
    if isfield(settings,'jackknife') & settings.jackknife == true
        Predict = false(Ndata,1);
        if strcmp(settings.method,'<=')
            Predict(data <= Result.Threshold) = true;
        elseif strcmp(settings.method,'<')
            Predict(data < Result.Threshold) = true;
        elseif strcmp(settings.method,'>=')
            Predict(data >= Result.Threshold) = true;
        else
            Predict(data > Result.Threshold) = true;
        end
        Rjack = lab_jackknife([outcome Predict],@ROCevaluate);
        Result.JACKKNIFE.Accuracy = Rjack.value(1,1);
        Result.JACKKNIFE.Sensitivity = Rjack.value(1,2);
        Result.JACKKNIFE.Specificity = Rjack.value(1,3);
        Result.JACKKNIFE.PPV = Rjack.value(1,4);
        Result.JACKKNIFE.NPV = Rjack.value(1,5);
        Result.JACKKNIFE.Precision = Rjack.value(1,6);
        Result.JACKKNIFE.F = Rjack.value(1,7);
        Result.JACKKNIFE.Gmean = Rjack.value(1,8);
        Result.JACKKNIFE.Accuracy_ci = [Rjack.value_low(1,1) Rjack.value_high(1,1)];
        Result.JACKKNIFE.Sensitivity_ci = [Rjack.value_low(1,2) Rjack.value_high(1,2)];
        Result.JACKKNIFE.Specificity_ci = [Rjack.value_low(1,3) Rjack.value_high(1,3)];
        Result.JACKKNIFE.PPV_ci = [Rjack.value_low(1,4) Rjack.value_high(1,4)];
        Result.JACKKNIFE.NPV_ci = [Rjack.value_low(1,5) Rjack.value_high(1,5)];
        Result.JACKKNIFE.Precision_ci = [Rjack.value_low(1,6) Rjack.value_high(1,6)];
        Result.JACKKNIFE.F_ci = [Rjack.value_low(1,7) Rjack.value_high(1,7)];
        Result.JACKKNIFE.Gmean_ci = [Rjack.value_low(1,8) Rjack.value_high(1,8)];
        clearvars Rjack
        if strcmp(settings.method,'<=') | strcmp(settings.method,'<')
            Rjack = lab_jackknife([outcome -data],@AUCevaluate);
        else
            Rjack = lab_jackknife([outcome data],@AUCevaluate);
        end
        Result.JACKKNIFE.AUC = Rjack.value(1,1);
        Result.JACKKNIFE.AUC_ci = [Rjack.value_low(1,1) Rjack.value_high(1,1)];
    end
    
    % Plot Result
    if ~isfield(settings,'filename') | ~isempty(settings.filename)
        fig1 = figure('NumberTitle','off','Name','ROC-Curve','MenuBar','none','Color',[1 1 1]);
        m1 = uimenu(fig1,'Label','File');
        uimenu(m1,'Label','Save Picture','Callback',@(~,~)lab_print_figure);
        uimenu(m1,'Label','Copy Picture','Callback','print -dmeta;');
        uimenu(m1,'Label','Close','Callback','close;');
        Ncolor = [1 0.9 0.9];
        if isfield(Result,'BOOTSTRAP') & isfield(Result.BOOTSTRAP,'Sensitivity')
            Rw = Result.BOOTSTRAP.Sensitivity(2) - Result.BOOTSTRAP.Sensitivity(1);
            Rh = Result.BOOTSTRAP.Specificity(2)  - Result.BOOTSTRAP.Specificity(1);
            rectangle('Position',[Result.BOOTSTRAP.Specificity(1) Result.BOOTSTRAP.Sensitivity(1) Rh Rw],'Facecolor',Ncolor,'EdgeColor','none');
            hold on
            Ncolor = [0.9 1 0.9];
        end
        if isfield(Result,'JACKKNIFE') & isfield(Result.JACKKNIFE,'Sensitivity_ci')
            Rw = Result.JACKKNIFE.Sensitivity_ci(2) - Result.JACKKNIFE.Sensitivity_ci(1);
            Rh = Result.JACKKNIFE.Specificity_ci(2)  - Result.JACKKNIFE.Specificity_ci(1);
            rectangle('Position',[Result.JACKKNIFE.Specificity_ci(1) Result.JACKKNIFE.Sensitivity_ci(1) Rh Rw],'Facecolor',Ncolor,'EdgeColor','none');
            hold on
        end
        plot(Result.Specificity,Result.Sensitivity)
        set(gca,'XDir','reverse','XTickLabel',num2str((0:10:100)'),'XLim',[0 1],'YTickLabel',num2str((0:10:100)'),'YLim',[0 1]);
        xlabel('Specificity');
        ylabel('Sensitivity');
        hold on
        scatter(Result.MaxYouden.Specificity,Result.MaxYouden.Sensitivity,'r');
        
        if isfield(settings,'filename')
            [~,ROC_filepath,~,ROC_fileS] = lab_filename(settings.filename);
            if isempty(ROC_filepath)
                ROC_filepath = pwd;
            end
            ROC_file = fullfile(ROC_filepath,[ROC_fileS '_ROC.tif']);
            lab_print_figure(ROC_file,fig1);
            
            xlsout = {'ROC-Curve','Max Youden-Index'};
            xlsout(2,1:2) = {'Threshold',Result.MaxYouden.Threshold};
            xlsout(3,1:2) = {'Sensitivity',Result.MaxYouden.Sensitivity};
            xlsout(4,1:2) = {'Specificity',Result.MaxYouden.Specificity};
            xlsout(5,1:2) = {'PPV',Result.MaxYouden.PPV};
            xlsout(6,1:2) = {'NPV',Result.MaxYouden.NPV};
            xlsout(7,1:2) = {'AUC',Result.AUC};
            if isfield(Result,'BOOTSTRAP') & isfield(Result.BOOTSTRAP,'Sensitivity')
                xlsout(1,end+1:end+2) = {'95% CI-Low','95% CI-High'};
                xlsout(3,end-1:end) = num2cell(Result.BOOTSTRAP.Sensitivity);
                xlsout(4,end-1:end) = num2cell(Result.BOOTSTRAP.Specificity);
                xlsout(5,end-1:end) = num2cell(Result.BOOTSTRAP.PPV);
                xlsout(6,end-1:end) = num2cell(Result.BOOTSTRAP.NPV);
                xlsout(7,end-1:end) = num2cell(Result.BOOTSTRAP.AUC);
            end
            if isfield(Result,'JACKKNIFE') & isfield(Result.JACKKNIFE,'Sensitivity_ci')
                xlsout(1,end+1:end+2) = {'Jackknife-Low','Jackknife-High'};
                xlsout(3,end-1:end) = num2cell(Result.JACKKNIFE.Sensitivity_ci);
                xlsout(4,end-1:end) = num2cell(Result.JACKKNIFE.Specificity_ci);
                xlsout(5,end-1:end) = num2cell(Result.JACKKNIFE.PPV_ci);
                xlsout(6,end-1:end) = num2cell(Result.JACKKNIFE.NPV_ci);
                xlsout(7,end-1:end) = num2cell(Result.JACKKNIFE.AUC_ci);
            end
            XLS_file = fullfile(ROC_filepath,[ROC_fileS '_ROC.xls']);
            lab_write_xls(XLS_file,xlsout);
        end
    end
end

function R = ROCevaluate(data)
    outcome = data(:,1);
    predict = data(:,2);
    
    idx = (outcome()==1);
    
    p = length(outcome(idx));
    n = length(outcome(~idx));
    N = p+n;
    
    tp = sum(outcome(idx)==predict(idx));
    tn = sum(outcome(~idx)==predict(~idx));
    fp = n-tn;
    fn = p-tp;
    
    Accuracy = (tp+tn)/N;
    Sensitivity = tp/p;
    Specificity = tn/n;
    PPV = tp / (tp + fp);
    NPV = tn / (tn + fn);
    Precision = tp/(tp+fp);
    F = 2*((Precision*Sensitivity)/(Precision + Sensitivity));
    Gmean = sqrt(Sensitivity*Specificity);
    
    R = [Accuracy Sensitivity Specificity PPV NPV Precision F Gmean];
end

function AUC = AUCevaluate(data)
    outcome = logical(data(:,1));
    predict = data(:,2);
    [~,~,~,AUC] = perfcurve(outcome,predict,true);
end