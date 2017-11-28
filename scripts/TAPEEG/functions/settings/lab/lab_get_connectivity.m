function [settings,skipprocessing] = lab_get_connectivity(settings,data,header,cfg)

skipprocessing = 0;

if length(setdiff({'PLI','dPLI','wPLI','PLT','PLV','wPLV','PTE','rPTE'},settings.measure)) < 8
    if ~isfield(settings,'PHASE') | isempty(settings.PHASE)
        settings = lab_get_PHASE(settings);
    end
else
    settings.PHASE = [];
end

if max(strcmp(settings.measure,'SL'))
    if ~isfield(settings,'SL') | isempty(settings.SL)
        settings = lab_get_SL(settings);
    end
else
    settings.SL = [];
end

if max(strcmp(settings.measure,'SLc'))
    if ~isfield(settings,'SLc') | isempty(settings.SLc)
        settings = lab_get_SLc(settings);
    end
else
    settings.SLc = [];
end

if max(strcmp(settings.measure,'EC'))
    if ~isfield(settings,'EC') | isempty(settings.EC)
        settings = lab_get_EC(settings);
    end
else
    settings.EC = [];
end

if max(strcmp(settings.measure,'DTF'))
    if ~isfield(settings,'DTF') | isempty(settings.DTF)
        settings = lab_get_MARKERS(settings,header,cfg);
        settings = lab_get_DTF(settings,data,header);
    end
else
    settings.DTF = [];
end

if max(strcmp(settings.measure,'TE'))
    if ~isfield(settings,'TE') | isempty(settings.TE)
        settings = lab_get_TE(settings);
    end
else
    settings.TE = [];
end