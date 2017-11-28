function [cfg,skipprocessing] =  lab_set_svm(cfg)

skipprocessing = 0;

if ~exist('cfg','var') | ~isfield(cfg,'SVM') | ~isfield(cfg.SVM,'kernel')
    cfg.SVM.kernel = 'gaussianKernel';
    cfg.SVM.regularization = 1;
    cfg.SVM.sigma = 0.1;
    cfg.SVM.autosigma = false;
    cfg.SVM.max_passes = 5;
    cfg.SVM.tol = 0.001;
    cfg.SVM.percenttest = 0.3;
    cfg.SVM.loops = 100;
end

Prompt = cell(0,2);
Formats = {};
Prompt(end+1,:) = {'Select kernel','kernel'};
Formats(end+1,1).type = 'list';
Formats(end,1).style = 'popupmenu';
Formats(end,1).format = 'input';
Formats(end,1).items = {'gaussianKernel','linearKernel'};

Prompt(end+1,:) = {'Automatic sigma (gaussian kernel)' 'autosigma'};
Formats(end+1,1).type = 'check';

[cfg.SVM,Cancelled] = inputsdlg(Prompt,'SVM',Formats,cfg.SVM);
if Cancelled == 1
    skipprocessing = 1;
    return
end

if strcmp(cfg.SVM.kernel,'linearKernel')
    cfg.SVM.autosigma = false;
end

Prompt = cell(0,2);
Formats = {};
if strcmp(cfg.SVM.kernel,'gaussianKernel') & cfg.SVM.autosigma == false
    Prompt(end+1,:) = {'SVM regularization', 'regularization'};
    Formats(end+1,1).type = 'edit';
    Formats(end,1).format = 'integer';
    Formats(end,1).limits = [0 99];
    Formats(end,1).size = 40;
    
    Prompt(end+1,:) = {'Sigma (gaussian kernel)', 'sigma'};
    Formats(end+1,1).type = 'edit';
    Formats(end,1).format = 'float';
    Formats(end,1).limits = [0 inf];
    Formats(end,1).size = 40;
elseif strcmp(cfg.SVM.kernel,'gaussianKernel') & cfg.SVM.autosigma == true
    Prompt(end+1,:) = {'Percent of trials for testing automatic sigma', 'percenttest'};
    Formats(end+1,1).type = 'edit';
    Formats(end,1).format = 'float';
    Formats(end,1).limits = [0 inf];
    Formats(end,1).size = 40;
    
    Prompt(end+1,:) = {'Number of loops for testing automatic sigma', 'loops'};
    Formats(end+1,1).type = 'edit';
    Formats(end,1).format = 'integer';
    Formats(end,1).limits = [0 99999];
    Formats(end,1).size = 40;
elseif strcmp(cfg.SVM.kernel,'linearKernel')
    Prompt(end+1,:) = {'SVM regularization', 'regularization'};
    Formats(end+1,1).type = 'edit';
    Formats(end,1).format = 'integer';
    Formats(end,1).limits = [0 99];
    Formats(end,1).size = 40;
        
    Prompt(end+1,:) = {'Number of iterations', 'max_passes'};
    Formats(end+1,1).type = 'edit';
    Formats(end,1).format = 'integer';
    Formats(end,1).limits = [0 99];
    Formats(end,1).size = 40;
    
    Prompt(end+1,:) = {'Tolerance value (floating point)', 'tol'};
    Formats(end+1,1).type = 'edit';
    Formats(end,1).format = 'float';
    Formats(end,1).limits = [0 inf];
    Formats(end,1).size = 40;
end

if ~isempty(Prompt)
    [cfg.SVM,Cancelled] = inputsdlg(Prompt,'SVM',Formats,cfg.SVM);
    if Cancelled == 1
        skipprocessing = 1;
        return
    end
end