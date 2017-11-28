function settings = lab_get_PHASE(settings)

if exist('settings','var') & isfield(settings,'measure') & ...
        length(setdiff({'PLI','dPLI','wPLI','PLT','PLV','wPLV','PTE','rPTE'},settings.measure)) == 8
    settings.PHASE = [];
    return
end

if ~exist('settings','var') | ~isfield(settings,'PHASE') | ~isfield(settings.PHASE,'phaseestimate')
    settings.PHASE.phaseestimate = 'hilbert';
    settings.PHASE.window = 4;
    settings.PHASE.usehann = true;
    settings.PHASE.storephase = false;
    settings.PHASE.maxerror = 1;
    settings.PHASE.storeangle = false;
    settings.PHASE.debiased = true;
    settings.PHASE.step = 4096;
    settings.PHASE.stepunit = 'TFs';
    settings.PHASE.dosinus = true;
    settings.PHASE.delta = 2;
    settings.PHASE.deltafreq = true;
    settings.PHASE.dodirect = false;
    settings.PHASE.donorm = true;
    settings.PHASE.downsample = 2;
end

Prompt = cell(0,2);
Formats = [];

Prompt{end+1,1} = 'Phase estimation:';
Formats(end+1,1).type = 'text';
Formats(end,1).span = [1 2];

Prompt(end+1,:) = {'Phase estimate','phaseestimate'};
Formats(end+1,1).type = 'list';
Formats(end,1).style = 'popupmenu';
Formats(end,1).format = 'input';
Formats(end,1).items = {'hilbert';'zerocross';'wavelet';'rihaczek';'realphase'};
Formats(end,1).callback = {@correct_phaseestimate,'@ALL','@ALL'};

Prompt(end+1,:) = {'Use 50% Hanning window','usehann'};
Formats(end+1,1).type = 'check';
Formats(end,1).format = 'integer';
Formats(end,1).callback = {@correct_hann,'@ALL','@ALL'};

Prompt(end+1,:) = {'Sliding window length (seconds / 0 = no window)','window'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'float';
Formats(end,1).limits = [0 inf];
Formats(end,1).size = 60;
Formats(end,1).span = [1 2];

Prompt(end+1,:) = {'Skip timeframes at start','skipstart'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'integer';
Formats(end,1).limits = [0 inf];
Formats(end,1).size = 60;

Prompt(end+1,:) = {'at end','skipend'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'integer';
Formats(end,1).limits = [0 inf];
Formats(end,1).size = 60;

Prompt(end+1,:) = {'Store phase','storephase'};
Formats(end+1,1).type = 'check';
Formats(end,1).span = [1 2];

% Prompt(end+1,:) = {'Round <10^-12 angles (deprecated)','flagold'};
% Formats(end+1,1).type = 'check';
% Formats(end,1).span = [1 2];

Prompt{end+1,1} = ' ';
Formats(end+1,1).type = 'text';
Formats(end,1).span = [1 2];

Prompt{end+1,1} = 'Connectivity calculation:';
Formats(end+1,1).type = 'text';
Formats(end,1).span = [1 2];

Prompt(end+1,:) = {'Stepsize (0 = no steps)','step'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'float';
Formats(end,1).limits = [0 9999];
Formats(end,1).size = 40;

Prompt(end+1,:) = {'','stepunit'};
Formats(end+1,1).type = 'list';
Formats(end,1).style = 'popupmenu';
Formats(end,1).format = 'input';
Formats(end,1).items = {'phases of max frequency';'TFs';'seconds'};

if ~isfield(settings,'measure') | any(strcmp(settings.measure,'PLI')) | any(strcmp(settings.measure,'dPLI'))
    Prompt{end+1,1} = ' ';
    Formats(end+1,1).type = 'text';
    Formats(end,1).span = [1 2];
    
    Prompt{end+1,1} = '(d)PLI settings:';
    Formats(end+1,1).type = 'text';
    Formats(end,1).span = [1 2];
    
    Prompt(end+1,:) = {'Use sinus function','dosinus'};
    Formats(end+1,1).type = 'check';
    Formats(end,1).span = [1 2];
end

if ~isfield(settings,'measure') | any(strcmp(settings.measure,'wPLI'))
    Prompt{end+1,1} = ' ';
    Formats(end+1,1).type = 'text';
    Formats(end,1).span = [1 2];
    
    Prompt{end+1,1} = 'wPLI settings:';
    Formats(end+1,1).type = 'text';
    Formats(end,1).span = [1 2];
    
    Prompt(end+1,:) = {'Debiased','debiased'};
    Formats(end+1,1).type = 'check';
    Formats(end,1).span = [1 2];
end

if ~isfield(settings,'measure') | any(strcmp(settings.measure,'wPLV'))
    Prompt{end+1,1} = ' ';
    Formats(end+1,1).type = 'text';
    Formats(end,1).span = [1 2];
    
    Prompt{end+1,1} = 'wPLV settings:';
    Formats(end+1,1).type = 'text';
    Formats(end,1).span = [1 2];
    
    Prompt(end+1,:) = {'Max PhaseLag VolumeConduction','maxerror'};
    Formats(end+1,1).type = 'edit';
    Formats(end,1).format = 'float';
    Formats(end,1).limits = [0 20];
    Formats(end,1).size = 30;
    
    Prompt(end+1,:) = {'Store phase shift','storeangle'};
    Formats(end+1,1).type = 'check';
end

if ~isfield(settings,'measure') | any(strcmp(settings.measure,'PTE')) | any(strcmp(settings.measure,'rPTE'))
    Prompt{end+1,1} = ' ';
    Formats(end+1,1).type = 'text';
    Formats(end,1).span = [1 2];
    
    Prompt{end+1,1} = '(r)PTE settings:';
    Formats(end+1,1).type = 'text';
    Formats(end,1).span = [1 2];
    
    Prompt(end+1,:) = {'Delta','delta'};
    Formats(end+1,1).type = 'edit';
    Formats(end,1).format = 'float';
    Formats(end,1).limits = [0 inf];
    Formats(end,1).size = 30;
    
    Prompt(end+1,:) = {'Factor of median frequency','deltafreq'};
    Formats(end+1,1).type = 'check';
    Formats(end,1).format = 'input';
    Formats(end,1).callback = {@delta_frequency,'@ALL','@ALL'};
    
    if any(strcmp(settings.measure,'rPTE'))
        Prompt(end+1,:) = {'Downsample (Nyquist rate / n, 0 = disabled)','downsample'};
        Formats(end+1,1).type = 'edit';
        Formats(end,1).format = 'float';
        Formats(end,1).limits = [0 inf];
        Formats(end,1).size = 30;
    end
    
    Prompt(end+1,:) = {'Make directional','dodirect'};
    Formats(end+1,1).type = 'check';
    Formats(end,1).format = 'input';
    Formats(end,1).callback = {@pte_direct,'@ALL','@ALL'};
    Formats(end,1).span = [1 2];
    
    Prompt(end+1,:) = {'Normalize','donorm'};
    Formats(end+1,1).type = 'check';
    Formats(end,1).format = 'input';
    Formats(end,1).callback = {@pte_norm,'@ALL','@ALL'};
    Formats(end,1).span = [1 2];
end

[settings.PHASE,Cancelled] = inputsdlg(Prompt,'PLI/PLV/PLT/PTE',Formats,settings.PHASE);
if isempty(settings.PHASE) | Cancelled == 1
    settings.PHASE = [];
    return
end

end

function settings = correct_phaseestimate(settings)
    if strcmp(settings.phaseestimate,'wavelet') | strcmp(settings.phaseestimate,'rihaczek')
        settings.usehann = false;
        settings.markerexclude = {};
        settings.storephase = false;
        settings.skipstart = [];
        settings.skipend = [];
    end
end

function settings = correct_hann(settings)
    if strcmp(settings.phaseestimate,'wavelet') | strcmp(settings.phaseestimate,'rihaczek')
        settings.usehann = false;
        return
    end
    if settings.usehann == false
        settings.usehann = true;
    else
        settings.usehann = false;
    end
end

function settings = delta_frequency(settings)
    if settings.deltafreq == false
        settings.deltafreq = true;
        settings.delta = 2;
    else
        settings.deltafreq = false;
        settings.delta = 40;
    end
end

function settings = pte_direct(settings)
    if settings.donorm == true
        settings.dodirect = false;
    elseif settings.dodirect == true
        settings.dodirect = false;
    else
        settings.dodirect = true;
    end
end

function settings = pte_norm(settings)
    if settings.donorm == false
        settings.donorm = true;
        settings.dodirect = false;
    else
        settings.donorm = false;
    end
end