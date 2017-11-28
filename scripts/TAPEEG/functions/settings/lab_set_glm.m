function [cfg,skipprocessing] = lab_set_glm(cfg,header,result)

skipprocessing = 0;

if ~exist('header','var')
    header.vars = {'var1','var2'};
    header.result = {'result'};
end
if ~exist('result','var')
    result = [];
end

cfg.GLM.folder = 'GLM';
cfg.GLM.select = [];
if isfield(header,'vars')
    cfg.GLM.outcome = length(header.vars);
else
    cfg.GLM.outcome = [];
end
cfg.GLM.method = 'fit';
cfg.GLM.bootstrap = 0;
cfg.GLM.model = 'constant';
try %#ok<TRYNC>
    if length(unique(result(:,1))) == 2
        cfg.GLM.distribution = 'binomial';
        cfg.GLM.link = 'logit';
    elseif islogical(result(:,1))
        cfg.GLM.distribution = 'binomial';
        cfg.GLM.link = 'logit';
    end
    clearvars tmp
end
cfg.GLM.ROC.method = '>=';
cfg.GLM.ROC.nbootstrap = 0;
cfg.GLM.ROC.jackknife = false;

Prompt = {};
Formats = {};

Prompt(end+1,:) = {'Output-folder', 'folder'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'text';
Formats(end,1).size = 100;
Formats(end,1).span = [1 2];

Prompt(end+1,:) = {'Select measures','select'};
Formats(end+1,1).type = 'list';
Formats(end,1).style = 'listbox';
Formats(end,1).items = [header.vars header.result];
Formats(end,1).limits = [0 4]; % multi-select
Formats(end,1).size = [140 300];
Formats(end,1).span = [7 2];

Formats(end+1,1).type = 'none';
Formats(end,1).span = [1 2];

Prompt(end+1,:) = {'Select outcome','outcome'};
Formats(end+1,1).type = 'list';
Formats(end,1).style = 'popuplist';
Formats(end,1).items = header.result;
Formats(end,1).callback = {@switchoutcome,{'distribution','link'}, ...
    'outcome','distribution','link',result};
Formats(end,1).span = [1 2];

Formats(end+1,1).type = 'none';
Formats(end,1).span = [1 2];

Prompt(end+1,:) = {'Select method','method'};
Formats(end+1,1).type = 'list';
Formats(end,1).style = 'popuplist';
Formats(end,1).format = 'input';
Formats(end,1).items = {'fit','stepwise'};
Formats(end,1).callback = {@set_method,'bootstrap','method','bootstrap'};

Prompt(end+1,:) = {'Bootstrap (0 = off)','bootstrap'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'integer';
Formats(end,1).limits = [0 999999];
Formats(end,1).size = 50;

Formats(end+1,1).type = 'none';
Formats(end,1).span = [1 2];

Prompt(end+1,:) = {'Select distribution','distribution'};
Formats(end+1,1).type = 'list';
Formats(end,1).style = 'popuplist';
Formats(end,1).format = 'input';
Formats(end,1).items = {'binomial','gamma','inverse gaussian','normal','poisson'};

Prompt(end+1,:) = {'Select link function','link'};
Formats(end+1,1).type = 'list';
Formats(end,1).style = 'popuplist';
Formats(end,1).format = 'input';
Formats(end,1).items = {'identity','log','logit','probit','comploglog','reciprocal','loglog'};

Formats(end+1,1).type = 'none';
Formats(end,1).span = [1 2];

Prompt(end+1,:) = {'Select model','model'};
Formats(end+1,1).type = 'list';
Formats(end,1).style = 'popuplist';
Formats(end,1).format = 'input';
Formats(end,1).items = {'constant','linear','interactions','purequadratic','quadratic','polyijk'};

Prompt(end+1,:) = {'Select criterion','criterion'};
Formats(end+1,1).type = 'list';
Formats(end,1).style = 'popuplist';
Formats(end,1).format = 'input';
Formats(end,1).items = {'sse','deviance','aic','bic','rsquared','adjrsquared'};

Formats(end+1,1).type = 'none';
Formats(end,1).span = [1 2];

Prompt(end+1,:) = {'ROC','ROC'};
Formats(end+1,1).type = 'check';
Formats(end,1).callback = {@lab_set_create_ROC,'ROC','ROC',result};

[cfg.GLM,Cancelled] = inputsdlg(Prompt,'Generalized linear Regression',Formats,cfg.GLM);
if Cancelled == 1 | isempty(cfg.GLM.select)
    cfg.GLM = [];
    skipprocessing = 1;
    return
end

end

function [distribution,link] = switchoutcome(outcome,distribution,link,result)
    try %#ok<TRYNC>
        tmp = result(:,outcome);
        if length(unique(tmp)) == 2
            distribution = 'binomial';
            link = 'logit';
        elseif islogical(tmp)
            distribution = 'binomial';
            link = 'logit';
        else
            distribution = 'normal';
            link = 'identity';
        end
        clearvars tmp
    end
end

function bootstrap = set_method(method,bootstrap)
   if ~strcmp(method,'fit')
       bootstrap = 0;
   end
end