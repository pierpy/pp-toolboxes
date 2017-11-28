% Calculate Support Vector Machine
%
% [Result,cfg] = lab_svm(cfg)
%
% cfg (optional)
%
% Written by F. Hatz 2014

function [Result,cfg] = lab_svm(cfg)

Result = [];
if ~exist('cfg','var')
    cfg = [];
end
rng('default');
rng('shuffle');

[data,header,results,~,cfg] = lab_read_statistics(cfg,1,0);
if isempty(data) | size(results,2) ~= 1
    disp('Abort only input-data is not supported')
    return
end

% generate Output-path and -file
[~,~,~,filenameS] = lab_filename(cfg.file);
warning off %#ok<WNOFF>
mkdir(fullfile(cfg.path,'SVM'));
warning on %#ok<WNON>
cfg.path = fullfile(cfg.path,'SVM');

classes = unique(results);
if length(classes) ~= 2
    disp('Abort only outcome with two classes is supported')
    return
elseif isnumeric(results)
    classes = sort(classes);
    results(results == classes(1)) = 0;
    results(results == classes(2)) = 1;
else
    return
end

[cfg,skipprocessing] =  lab_set_svm(cfg);
if skipprocessing == 1
    return
else
    cfg.SVM = cfg.SVM;
end

ntrials = size(data,1);
for nResult = 1:size(results,2)
    result = results(:,nResult);
    
    filenameout = fullfile(cfg.path,[filenameS header.result{nResult}]);
    filenameS2 = [filenameS header.result{nResult}];
    
    % find sigma & regularization if necessary
    if cfg.SVM.autosigma == true 
        [cfg.SVM.regularization,cfg.SVM.sigma] = lab_example3parameters(data,result,cfg);
    end
    C = cfg.SVM.regularization;
    sigma = cfg.SVM.sigma;
    
    for ntest = 1:ntrials
        testvector = setdiff(1:ntrials,ntest);
        X_trn = data(testvector,:);
        Y_trn = result(testvector,1);
        X_tst = data(ntest,:);
        Y_tst = result(ntest,1);
        if strcmp(cfg.SVM.kernel,'linearKernel')
            modeltst = svmTrain(data,result,C,@linearKernel,cfg.SVM.tol,cfg.SVM.max_passes);
        else
            modeltst = svmTrain(data,result,C,@(x1, x2) gaussianKernel(x1, x2,sigma));
        end
        Predict_tst(ntest) = svmPredict(modeltst,X_tst);
    end
    
    % calculate full model
    if strcmp(cfg.SVM.kernel,'linearKernel')
        disp('Training Linear SVM ...')
        model = svmTrain(data,result,C,@linearKernel,cfg.SVM.tol,cfg.SVM.max_passes);
    else
        disp('Training SVM with RBF Kernel ...')
        model= svmTrain(data,result,C,@(x1, x2) gaussianKernel(x1, x2,sigma));
    end
    
    % Calculate Specificity and Sensitifity
    Predict_full = svmPredict(model,data);
    Sens_full = sum(Predict_full(result==1)) / sum(result);
    Spez_full = sum(Predict_full(result==0) == 0) / sum(result==0);
    Sens_tst = sum(Predict_tst(result==1)) / sum(result);
    Spez_tst = sum(Predict_tst(result==0) == 0) / sum(result==0);
    result = logical(result);
    [X,Y,~,AUC,OPTROCPT] = perfcurve(result,Predict_tst,true);
    Result.ROC.AUC = AUC;
    Result.ROC.X = X;
    Result.ROC.Y = Y;
    Result.ROC.sensitivity = OPTROCPT(2);
    Result.ROC.specificity = 1 - OPTROCPT(1);
    
    % write results
    if exist('filenameout','var') & ~isempty(filenameout)
        xlsout = {'Specificity_Full',Spez_full; ...
            'Sensitivity_Full',Sens_full; ...
            'Specificity_Test',Spez_tst; ...
            'Sensitivity_Test',Sens_tst; ...
            'B',model.b};
        for i = 1:size(data,2)
            xlsout(end+1,:) = {['W_' header.vars{i}],model.w(i)}; %#ok<AGROW>
        end
        lab_write_xls([filenameout '_SVM_' cfg.SVM.kernel '.xls'],xlsout);
        if isfield(model,'alphas')
            lab_write_xls([filenameout '_SVM_' cfg.SVM.kernel '_alphas.xls'],model.alphas);
        end
    end
    
    % plot if only 2 variables
    if size(X_trn,2) == 2
        fig1 = figure('Name','Decision boundary - train','MenuBar','none','Color',[1 1 1],'NumberTitle','off');
        m1 = uimenu(fig1,'Label','File');
        uimenu(m1,'Label','Save Picture','Callback',@(~,~)lab_print_figure);
        uimenu(m1,'Label','Copy Picture','Callback','print -dmeta;');
        uimenu(m1,'Label','Close','Callback','close;');
        visualizeBoundary(X_trn,Y_trn,model);
        title('Decision boundary - train');
        if exist('filenameout','var') & ~isempty(filenameout)
            lab_print_figure([filenameout '_' cfg.SVM.kernel '_DecisionBoundaryTrain.tif'],fig1);
        end
        
        fig2 = figure('Name','Decision boundary - test','MenuBar','none','Color',[1 1 1],'NumberTitle','off');
        m2 = uimenu(fig2,'Label','File');
        uimenu(m2,'Label','Save Picture','Callback',@(~,~)lab_print_figure);
        uimenu(m2,'Label','Copy Picture','Callback','print -dmeta;');
        uimenu(m2,'Label','Close','Callback','close;');
        visualizeBoundary(X_tst,Y_tst,model);
        title('Decision boundary - test');
        if exist('filenameout','var') & ~isempty(filenameout)
            lab_print_figure([filenameout '_' cfg.SVM.kernel '_DecisionBoundaryTest.tif'],fig2);
        end
    end
    
    % plot weights
    [~,tmp] = sort(model.w,'descend');
    weights = model.w(tmp);
    vars = header.vars(tmp);
    for i = 1:length(vars)
        vars{i} = [vars{i} ' '];
    end
    fig3 = figure('Name','Weights','MenuBar','none','Color',[1 1 1],'NumberTitle','off');
    m3 = uimenu(fig3,'Label','File');
    uimenu(m3,'Label','Save Picture','Callback',@(~,~)lab_print_figure);
    uimenu(m3,'Label','Copy Picture','Callback','print -dmeta;');
    uimenu(m3,'Label','Close','Callback','close;');
    plot(weights,'o','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',7);
    set(gca,'XTick',1:length(vars));
    set(gca,'XTickLabel',vars,'FontName','Times','fontsize',9);
    rotateXLabels(gca,90);
    title('Weights');
    if exist('filenameout','var') & ~isempty(filenameout)
        lab_print_figure([filenameout '_' cfg.SVM.kernel '_Weights.tif'],fig3);
    end
    
    % plot ROC
    fig4 = figure('Name','ROC','MenuBar','none','Color',[1 1 1]);
    m4 = uimenu(fig4,'Label','File');
    uimenu(m4,'Label','Save Picture','Callback',@(~,~)lab_print_figure);
    uimenu(m4,'Label','Copy Picture','Callback','print -dmeta;');
    uimenu(m4,'Label','Close','Callback','close;');
    plot(X,Y)
    xlabel('False positive rate'); ylabel('True positive rate')
    title(['ROC (AUC: ' num2str(AUC) ')']);
    hold on
    scatter(OPTROCPT(1),OPTROCPT(2),'r');
    if exist('filenameout','var') & ~isempty(filenameout)
        lab_print_figure([filenameout '_ROC.tif'],fig4);
    end
    
    % write verbose
    if exist('filenameout','var') & ~isempty(filenameout)
        fid=fopen([filenameout '_SVM.vrb'],'w');
        fprintf(fid,'SVM analysis\n');
        fprintf(fid,['File: ' filenameS2 '\n']);
        fprintf(fid,'\n');
        fprintf(fid,['Kernel: ' cfg.SVM.kernel '\n']);
        fprintf(fid,'\n');
        fprintf(fid,['Regularization: ' C '\n']);
        if strcmp(cfg.SVM.kernel,'gaussianKernel') & cfg.SVM.autosigma == false
            fprintf(fid,['Sigma: ' num2str(cfg.SVM.sigma) '\n']);
        elseif strcmp(cfg.SVM.kernel,'gaussianKernel') & cfg.SVM.autosigma == true
            fprintf(fid,['Sigma: ' num2str(sigma) '\n']);
        end
        fprintf(fid,['Number of iterations: ' num2str(cfg.SVM.max_passes) '\n']);
        fprintf(fid,['Tolerance value (floating point): ' num2str(cfg.SVM.tol) '\n']);
        fprintf(fid,['Percent trials for testing automatic sigma: ' num2str(cfg.SVM.percenttest) '\n']);
        fclose(fid);
    end
    
    Result.model = model;
end