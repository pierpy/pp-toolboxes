function [cfg,skipprocessing] = lab_set_randomforest(cfg,ntrials,nvars,nclass)

skipprocessing = 0;
if ~exist('ntrials','var')
    ntrials = [];
end
if ~exist('nvars','var')
    nvars = [];
end
if ~exist('nclass','var')
    nclass = [];
end

if ~exist('cfg','var') | ~isfield(cfg,'RANDFOREST') | ~isfield(cfg.RANDFOREST,'mode')
    if nclass > 5
        cfg.RANDFOREST.mode = 'regression';
    else
        cfg.RANDFOREST.mode = 'classification';
    end
end
Prompt = cell(0,2);
Formats = {};
Prompt(end+1,:) = {'Select mode','mode'};
Formats(end+1,1).type = 'list';
Formats(end,1).style = 'popupmenu';
Formats(end,1).format = 'input';
Formats(end,1).items = {'regression','classification'};

[cfg.RANDFOREST,Cancelled] = inputsdlg(Prompt,'Random Forest',Formats,cfg.RANDFOREST);
if Cancelled == 1
    skipprocessing = 1;
    return
end

if ~isfield(cfg.RANDFOREST,'ntree')
    cfg.RANDFOREST.nloops = 100;
    cfg.RANDFOREST.ntree = 1000;
    cfg.RANDFOREST.mtry = floor(sqrt(nvars));
    cfg.RANDFOREST.sampsize = ntrials;
    cfg.RANDFOREST.replace = true;
    if strcmp(cfg.RANDFOREST.mode,'classification')
        cfg.RANDFOREST.nodesize = 1;
    else
        cfg.RANDFOREST.nodesize = 5;
    end
    cfg.RANDFOREST.percenttest = 0;
    
    cfg.RANDFOREST.importance = false;
    cfg.RANDFOREST.localImp = false;
    cfg.RANDFOREST.proximity = false;
    cfg.RANDFOREST.oob_prox = true;
    if strcmp(cfg.RANDFOREST.mode,'regression')
        cfg.RANDFOREST.nPerm = 0;
        cfg.RANDFOREST.corr_bias = false;
    else
        cfg.RANDFOREST.classwt = ones(1,nclass);
        if ~isempty(nclass) & nclass > 0
            cfg.RANDFOREST.cutoff = ones(1,nclass)*(1/nclass);
        else
            cfg.RANDFOREST.cutoff = [];
        end
    end
end

Prompt = cell(0,2);
Formats = {};

Prompt(end+1,:) = {'Number of loops', 'nloops'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'integer';
Formats(end,1).limits = [0 inf];
Formats(end,1).size = 40;

Prompt(end+1,:) = {'Number trees', 'ntree'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'integer';
Formats(end,1).limits = [0 99999]; % 9-digits (positive #)
Formats(end,1).size = 60;

Prompt(end+1,:) = {'Mtry (default = sqrt(nvars))', 'mtry'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'integer';
Formats(end,1).limits = [0 99]; % 9-digits (positive #)
Formats(end,1).size = 40;

Prompt(end+1,:) = {'Size of sample to draw', 'sampsize'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'integer';
Formats(end,1).limits = [0 9999999]; % 9-digits (positive #)
Formats(end,1).size = 60;

Prompt(end+1,:) = {'Replacement' 'replace'};
Formats(end+1,1).type = 'check';
Formats(end,1).callback = {@setsampsize,'sampsize','replace',ntrials};

Prompt(end+1,:) = {'Minimum size of terminal nodes', 'nodesize'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'integer';
Formats(end,1).limits = [0 99]; % 9-digits (positive #)
Formats(end,1).size = 40;

Prompt{end+1,1} = ' ';
Formats(end+1,1).type = 'text';

Prompt(end+1,:) = {'Percent of trials for testing model (0 = all)', 'percenttest'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'float';
Formats(end,1).limits = [0 inf];
Formats(end,1).size = 40;

Prompt{end+1,1} = ' ';
Formats(end+1,1).type = 'text';

Prompt(end+1,:) = {'Importance of predictors' 'importance'};
Formats(end+1,1).type = 'check';

Prompt(end+1,:) = {'Casewise importance measure' 'localImp'};
Formats(end+1,1).type = 'check';

Prompt(end+1,:) = {'Proximity measure among the rows' 'proximity'};
Formats(end+1,1).type = 'check';

Prompt(end+1,:) = {'Proximity only on out-of-bag' 'oob_prox'};
Formats(end+1,1).type = 'check';

if isfield(cfg.RANDFOREST,'nPerm')
    Prompt(end+1,:) = {'Number of times the OOB data are permuted per tree', 'nPerm'};
    Formats(end+1,1).type = 'edit';
    Formats(end,1).format = 'integer';
    Formats(end,1).limits = [0 99]; % 9-digits (positive #)
    Formats(end,1).size = 40;
    
    Prompt(end+1,:) = {'Bias correction' 'corr_bias'};
    Formats(end+1,1).type = 'check';
else
    Prompt(end+1,:) = {'Weight of outcome classes','classwt'};
    Formats(end+1,1).type = 'edit';
    Formats(end,1).format = 'vector';
    Formats(end,1).limits = [-inf inf]; % if not [0 1] data must be numeric otherwise string
    
    Prompt(end+1,:) = {'Cutoff outcome classes','cutoff'};
    Formats(end+1,1).type = 'edit';
    Formats(end,1).format = 'vector';
    Formats(end,1).limits = [-inf inf]; % if not [0 1] data must be numeric otherwise string
    Formats(end,1).size = 300;
end

[cfg.RANDFOREST,Cancelled] = inputsdlg(Prompt,'Random Forest',Formats,cfg.RANDFOREST);
if Cancelled == 1
    skipprocessing = 1;
    return
end

end

function sampsize = setsampsize(replace,ntrials)
  if replace == true
      sampsize = ntrials;
  else
      sampsize = ceil(0.632*ntrials);
  end
end