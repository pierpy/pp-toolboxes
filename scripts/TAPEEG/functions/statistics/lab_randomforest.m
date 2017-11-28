% Calculate Random Forest Analysis
%
% [Result,cfg] = lab_randomforest(cfg)
%
% cfg (optional)
%
% Written by F. Hatz 2014

function [Result,cfg] = lab_randomforest(cfg)

Result = [];
if ~exist('cfg','var')
    cfg = [];
end
rng('default');
rng('shuffle');

[data,header,resultsAll,~,cfg] = lab_read_statistics(cfg,1,0);
if isempty(data)
    return
end

% generate Output-path and -file
[~,~,~,filenameS] = lab_filename(cfg.file);
warning off %#ok<WNOFF>
mkdir(fullfile(cfg.path,'RandomForest'));
warning on %#ok<WNON>
cfg.path = fullfile(cfg.path,'RandomForest');

% do settings
[ntrials,nvars] =size(data);
if size(resultsAll,2) == 1
    nclass = length(unique(resultsAll));
else
    nclass = [];
end
[cfg,skipprocessing] = lab_set_randomforest(cfg,ntrials,nvars,nclass);
pause(0.2);
if skipprocessing == 1;
    return
else
    settings = cfg.RANDFOREST;
end

for nResult = 1:size(resultsAll,2)
    filenameout = fullfile(cfg.path,[filenameS header.result{nResult}]);
    filenameS2 = [filenameS header.result{nResult}];
    
    % split in samples for model calculation and testing
    randvector = randperm(ntrials);
    cut = ceil(ntrials*(1-settings.percenttest));
    X_trn = data(randvector(1:cut),:);
    Y_trn = resultsAll(randvector(1:cut),nResult);
    if cut < size(data,1)
        X_tst = data(randvector(cut+1:end),:);
        Y_tst = resultsAll(randvector(cut+1:end),nResult);
    else
        X_tst = [];
        Y_tst = [];
    end
    settings.sampsize = floor(settings.sampsize * (size(X_trn,1)/size(data,1)));
    
    if isempty(settings.nloops) | settings.nloops < 1
        settings.nloops = 1;
    end
    
    importance = [];
    importanceSD = [];
    localImp = [];
    mse = [];
    for L = 1:settings.nloops
        % do calculation
        if strcmp(settings.mode,'regression')
            Model = regRF_train(X_trn,Y_trn,settings.ntree,settings.mtry,settings);
            if ~isempty(X_tst)
                Y_hat(:,L) = regRF_predict(X_tst,Model);
            end
        else
            Model = classRF_train(X_trn,Y_trn,settings.ntree,settings.mtry,settings);
            settings.predict_all = 1;
            if ~isempty(X_tst)
                [Y_hat(:,L),votes,prediction_pre_tree] = classRF_predict(X_tst,Model,settings);
            end
        end
        importance(:,:,L) = Model.importance;
        importanceSD(:,:,L) = Model.importanceSD;
        localImp(:,:,L) = Model.localImp;
        if isfield(Model,'mse')
            mse(:,:,L) = Model.mse;
        end
    end
    Model.importance = importance;
    Model.importanceSD = importanceSD;
    Model.localImp = localImp;
    Model.mse = mse;
    clearvars importance importanceSD localImp mse
    
    if ~isempty(X_tst)
        Y_tst = repmat(Y_tst,settings.nloops,1);
        Y_hat = reshape(Y_hat,size(Y_hat,1)*size(Y_hat,2),1);
    else
        alltrials = 1:size(data,1);
        settings.sampsize = settings.sampsize - 1;
        for L = 1:size(data,1)
            X_trn2 = X_trn(setdiff(alltrials,L),:);
            Y_trn2 = Y_trn(setdiff(alltrials,L),:);
            if strcmp(settings.mode,'regression')
                model2 = regRF_train(X_trn2,Y_trn2,settings.ntree,settings.mtry,settings);
                Y_hat(L,1) = regRF_predict(X_trn(L,:),model2);
            else
                model2 = classRF_train(X_trn2,Y_trn2,settings.ntree,settings.mtry,settings);
                settings.predict_all = 1;
                Y_hat(L,1) = classRF_predict(X_trn(L,:),model2,settings);
            end
        end
        Y_tst = Y_trn;
        settings.sampsize = settings.sampsize + 1;
    end
    
    if length(unique(Y_tst)) == 2
        [X,Y,~,AUC,OPTROCPT] = perfcurve(Y_tst,Y_hat,max(Y_tst));
        Result.AUC = AUC;
        Result.sensitivity = OPTROCPT(2);
        Result.specificity = 1 - OPTROCPT(1);
    end
    % collect results
    if strcmp(settings.mode,'regression')
        if ~isempty(Y_hat)
            Result.MSErate = sum((Y_hat-Y_tst).^2);
        else
            Result.MSErate= [];
        end
        Model.result = {'MSErate',Result.MSErate; ...
            'NumberTrees',Model.ntree;'mtry',Model.mtry; ...
            'Loops',settings.nloops};
    else
        if ~isempty(Y_hat)
            Result.error_rate = length(find(Y_hat~=Y_tst))/length(Y_tst);
        else
            Result.error_rate = [];
        end
        Model.result = {'Error rate',Result.error_rate; ...
            'NumberTrees',Model.ntree;'mtry',Model.mtry; ...
            'Loops',settings.nloops};
        if isfield(Result,'AUC')
            Model.result = cat(1,Model.result,{'','';'ROC','';'AUC',AUC; ...
                'Senstivity',OPTROCPT(2);'Specificity',1-OPTROCPT(1)});
        end
    end
    
    % write results to excel
    variables = {'result','importance','importanceSD','localImp','ntree','mtry',...
        'votes','oob_times','proximity','errtr','mse','rsq'};
    for i = 1:length(variables)
        if isfield(Model,variables{i})
            eval(['Result.' variables{i} '=Model.' variables{i} ';']);
            eval(['xlsout=Model.' variables{i} ';']);
            if size(xlsout,3) > 1 & isnumeric(xlsout)
                xlsout = mean(xlsout,3);
            end
            if exist('filenameout','var') & ~isempty(filenameout) & length(xlsout) > 1
                if size(xlsout,2) > 255
                    filenameout2 = [filenameout '_' variables{i} '.xlsx'];
                else
                    filenameout2 = [filenameout '_' variables{i} '.xls'];
                end
                lab_write_xls(filenameout2,xlsout);
                if strcmp(variables{i},'importance')
                    xlsout = {};
                    Dsort = Model.importance;
                    if size(Dsort,3) > 1
                        Dsort = mean(Dsort,3);
                    end
                    [~,Dsort] = sort(Dsort,'descend');
                    if isfield(cfg,'clustervars2') & cfg.clustervars2 > 1
                        xlsout{1,1} = ['C' num2str(cfg.clustervars) ' R' num2str(cfg.numresults) ' V' num2str(cfg.clustervars2)];
                    elseif isfield(cfg,'clustervars') & isfield(cfg,'numresults')
                        xlsout{1,1} = ['C' num2str(cfg.clustervars) ' R' num2str(cfg.numresults)];
                    else
                        xlsout{1,1} = '';
                    end
                    xlsout = [xlsout;header.subjects]; %#ok<AGROW>
                    datatmp = cat(1,header.vars,num2cell(data));
                    datatmp = datatmp(:,Dsort);
                    datatmp = [datatmp cat(1,header.result,num2cell(resultsAll))]; %#ok<AGROW>
                    xlsout = [xlsout datatmp]; %#ok<AGROW>
                    lab_write_xls([filenameout '_SortedByImportance.xlsx'],xlsout);
                end
            end
        end
    end
    
    % Plot results
    for i = 1:length(header.vars)
        varsplot{i} = regexprep(header.vars{i},'_',' '); %#ok<AGROW>
    end
    [~,sorttmp] = sort(median(Model.importance(:,end,:),3),'descend');
    fig1 = figure('Name','Importance Gini','MenuBar','none','Color',[1 1 1],'NumberTitle','off');
    m1 = uimenu(fig1,'Label','File');
    uimenu(m1,'Label','Save Picture','Callback',@(~,~)lab_print_figure);
    uimenu(m1,'Label','Copy Picture','Callback','print -dmeta;');
    uimenu(m1,'Label','Close','Callback','close;');
    boxplot(permute(Model.importance(sorttmp,end,:),[1 3 2])')
    set(gca,'XTick',1:length(Model.importance(:,end,1)),'XTickLabel',varsplot(sorttmp));
    rotateXLabels(gca,50);
    set(gca,'XTickLabel',repmat({' '},1,length(varsplot(sorttmp))));
    xlabel('feature');
    ylabel('magnitude');
    title('Mean decrease in Gini coefficient');
    if exist('filenameout','var') & ~isempty(filenameout)
        lab_print_figure([filenameout '_Importance.tif'],fig1);
    end
    if size(Model.importance,2) > 1
        [~,sorttmp] = sort(median(Model.importance(:,end-1,:),3),'descend');
        fig2 = figure('Name','Importance Accuracy','MenuBar','none','Color',[1 1 1],'NumberTitle','off');
        m2 = uimenu(fig2,'Label','File');
        uimenu(m2,'Label','Save Picture','Callback',@(~,~)lab_print_figure);
        uimenu(m2,'Label','Copy Picture','Callback','print -dmeta;');
        uimenu(m2,'Label','Close','Callback','close;');
        boxplot(permute(Model.importance(sorttmp,end-1,:),[1 3 2])')
        xlabel('feature');
        ylabel('magnitude');
        set(gca,'XTick',1:length(Model.importance(:,end,1)),'XTickLabel',varsplot(sorttmp));
        rotateXLabels(gca,50);
        title('Mean decrease in Accuracy');
        if exist('filenameout','var') & ~isempty(filenameout)
            lab_print_figure([filenameout '_ImportanceAccuracy.tif'],fig2);
        end
    end
    %importanceSD = The ?standard errors? of the permutation-based importance measure. For classification,
    %           a D by nclass + 1 matrix corresponding to the first nclass + 1
    %           columns of the importance matrix. For regression, a length p vector.
    if length(Model.importanceSD(:,:,1)) > 1
        fig3 = figure('Name','Importance SD','MenuBar','none','Color',[1 1 1],'NumberTitle','off');
        m3 = uimenu(fig3,'Label','File');
        uimenu(m3,'Label','Save Picture','Callback',@(~,~)lab_print_figure);
        uimenu(m3,'Label','Copy Picture','Callback','print -dmeta;');
        uimenu(m3,'Label','Close','Callback','close;');
        boxplot(permute(Model.importanceSD(sorttmp,end,:),[1 3 2])')
        xlabel('feature');
        ylabel('magnitude');
        set(gca,'XTick',1:length(Model.importance(:,end,1)),'XTickLabel',varsplot(sorttmp));
        rotateXLabels(gca,50);
        title('Std. errors of importance measure');
        if exist('filenameout','var') & ~isempty(filenameout)
            lab_print_figure([filenameout '_ImportanceSD.tif'],fig3);
        end
    end
    if isfield(Model,'localImp') & length(Model.localImp(:,1,1)) > 1
        fig4 = figure('Name','Local importance','MenuBar','none','Color',[1 1 1],'NumberTitle','off');
        m4 = uimenu(fig4,'Label','File');
        uimenu(m4,'Label','Save Picture','Callback',@(~,~)lab_print_figure);
        uimenu(m4,'Label','Copy Picture','Callback','print -dmeta;');
        uimenu(m4,'Label','Close','Callback','close;');
        plot(mean(Model.localImp,3)');
        legend(varsplot,'Location','EastOutside');
        if exist('filenameout','var') & ~isempty(filenameout)
            lab_print_figure([filenameout '_LocalImportance.tif'],fig4);
        end
    end
    if isfield(Model,'mse') & length(Model.mse) > 1
        fig5 = figure('Name','OOB error rate','MenuBar','none','Color',[1 1 1],'NumberTitle','off');
        m5 = uimenu(fig5,'Label','File');
        uimenu(m5,'Label','Save Picture','Callback',@(~,~)lab_print_figure);
        uimenu(m5,'Label','Copy Picture','Callback','print -dmeta;');
        uimenu(m5,'Label','Close','Callback','close;');
        plot(mean(Model.mse,3)); title('OOB MSE error rate');  xlabel('iteration (# trees)'); ylabel('OOB error rate');
        if exist('filenameout','var') & ~isempty(filenameout)
            lab_print_figure([filenameout '_OOBerror.tif'],fig5);
        end
    end
    
    if isfield(Result,'AUC')
        fig6 = figure('Name','ROC','MenuBar','none','Color',[1 1 1],'NumberTitle','off');
        m6 = uimenu(fig6,'Label','File');
        uimenu(m6,'Label','Save Picture','Callback',@(~,~)lab_print_figure);
        uimenu(m6,'Label','Copy Picture','Callback','print -dmeta;');
        uimenu(m6,'Label','Close','Callback','close;');
        plot(X,Y);
        xlabel('False positive rate');
        ylabel('True positive rate');
        title(['ROC for classification by logistic regression (AUC: ' num2str(Result.AUC) ')']);
        hold on
        scatter(OPTROCPT(1),OPTROCPT(2),'r');
        lab_print_figure([filenameout '_ROC.tif'],fig6);
    elseif ~isempty(Y_hat) & isfield(Result,'MSErate') & ~isempty(Result.MSErate)
        fig6 = figure('Name','Predicted-True Values','MenuBar','none','Color',[1 1 1],'NumberTitle','off');
        m6 = uimenu(fig6,'Label','File');
        uimenu(m6,'Label','Save Picture','Callback',@(~,~)lab_print_figure);
        uimenu(m6,'Label','Copy Picture','Callback','print -dmeta;');
        uimenu(m6,'Label','Close','Callback','close;');
        plot(Y_tst,Y_hat);
        xlabel('True values');
        ylabel('Predicted values'); 
        title(['Predicted versus true outcome values (MSE rate: ' num2str(Result.MSErate) ')']);
        lab_print_figure([filenameout '_MSE.tif'],fig6);
    end 
    
    % Save Matlab-Container with result
    save([filenameout '_Model.mat'],'Model','Result');
    
    % write verbose file
    if exist('filenameout','var') & ~isempty(filenameout)
        fid=fopen([filenameout '_randomforest.vrb'],'w');
        fprintf(fid,'Random forest analysis\n');
        fprintf(fid,['File: ' filenameS2 '\n']);
        fprintf(fid,'\n');
        fprintf(fid,['Method: ' cfg.RANDFOREST.mode '\n']);
        fprintf(fid,'\n');
        fprintf(fid,['Number of loops: ' num2str(cfg.RANDFOREST.nloops) '\n']);
        fprintf(fid,['Number of trees: ' num2str(cfg.RANDFOREST.ntree) '\n']);
        fprintf(fid,['Mtry: ' num2str(cfg.RANDFOREST.mtry) '\n']);
        fprintf(fid,['Sample size: ' num2str(cfg.RANDFOREST.sampsize) '\n']);
        fprintf(fid,['Replace: ' num2str(cfg.RANDFOREST.replace) '\n']);
        fprintf(fid,['Node size: ' num2str(cfg.RANDFOREST.nodesize) '\n']);
        fprintf(fid,['Importance: ' num2str(cfg.RANDFOREST.importance) '\n']);
        fprintf(fid,['Local importance: ' num2str(cfg.RANDFOREST.localImp) '\n']);
        fprintf(fid,['Proximity: ' num2str(cfg.RANDFOREST.proximity) '\n']);
        fprintf(fid,['out-of-box proximity: ' num2str(cfg.RANDFOREST.oob_prox) '\n']);
        if isfield(cfg.RANDFOREST,'classwt')
            fprintf(fid,['Classes weights: ' num2str(cfg.RANDFOREST.classwt) '\n']);
            fprintf(fid,['Classes cutoff: ' num2str(cfg.RANDFOREST.cutoff) '\n']);
        end
        if isfield(cfg.RANDFOREST,'nPerm')
            fprintf(fid,['Number of permutations: ' num2str(cfg.RANDFOREST.nPerm) '\n']);
            fprintf(fid,['Bias correction: ' num2str(cfg.RANDFOREST.corr_bias) '\n']);
        end
        fprintf(fid,['Percent trials for testing: ' num2str(cfg.RANDFOREST.percenttest) '\n']);
        fprintf(fid,'\n');
        fclose(fid);
    end
end

return

