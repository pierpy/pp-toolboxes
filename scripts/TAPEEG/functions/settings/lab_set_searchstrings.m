function [cfg,skipprocessing] = lab_set_searchstrings(cfg,dofolder,doextra,Strings)

skipprocessing = 0;

if ~exist('Strings','var') | ~iscell(Strings)
    Strings = {'strings','Cartool (.sef)','ICA','Filtered','Results','Cartool (.ris)', ...
        'Cartool (.eph)','EDF (.edf)','EGI (.mff)','EGI (.raw)','EGI (.mat)','Text (.txt)', ...
        'Elekta (.fif)','MRI (.hdr)','Matrix (.txt)','Matrix (.mat)',};
end
if ~exist('doextra','var')
    doextra = 0;
end
if ~exist('dofolder','var')
    dofolder = 0;
end

if ~exist('cfg','var') | ~isfield(cfg,'SEARCH')
    cfg.SEARCH = [];
end
if isfield(cfg.SEARCH,'searchstring')
    if ischar(cfg.SEARCH.searchstring)
        tmp = cfg.SEARCH.searchstring;
        cfg.SEARCH.searchstring = [];
        cfg.SEARCH.searchstring{1} = tmp;
        clearvars tmp
    elseif ~iscell(cfg.SEARCH.searchstring)
         cfg.SEARCH.searchstring{1} = '';
    end
    cfg.SEARCH = set_searchstring(cfg.SEARCH);
end
if isfield(cfg.SEARCH,'excludestring')
    if ischar(cfg.SEARCH.excludestring)
        tmp = cfg.SEARCH.excludestring;
        cfg.SEARCH.excludestring = [];
        cfg.SEARCH.excludestring{1} = tmp;
        clearvars tmp
    elseif ~iscell(cfg.SEARCH.excludestring)
         cfg.SEARCH.excludestring{1} = '';
    end
end
if ~isfield(cfg.SEARCH,'searchmode')
    cfg.SEARCH.searchmode = Strings{2};
end
if ~isfield(cfg.SEARCH,'searchstring') | isempty(cfg.SEARCH.searchstring)
    cfg.SEARCH = set_searchmode(cfg.SEARCH);
end
if ~isfield(cfg.SEARCH,'excludeprocessed')
    cfg.SEARCH.excludeprocessed = false;
end
if ~isfield(cfg.SEARCH,'doICAresult')
    cfg.SEARCH.doICAresult = false;
end

Formats = {};
Prompt = cell(0,2);

if dofolder == 1
    if ~isfield(cfg.SEARCH,'searchfolder') | ~ischar(cfg.SEARCH.searchfolder)
        cfg.SEARCH.searchfolder = '';
    end
    Prompt(end+1,:) = {'Search Folder','searchfolder'};
    Formats(end+1,1).type = 'edit';
    Formats(end,1).format = 'dir';
end

Prompt(end+1,:) = {'Search ','searchmode'};
Formats(end+1,1).type = 'list';
Formats(end,1).style = 'popupmenu';
Formats(end,1).format = 'input';
Formats(end,1).items = Strings;
Formats(end,1).callback = {@set_searchmode,'@ALL','@ALL'};

Prompt(end+1,:) = {'','searchstring'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'vector';
Formats(end,1).size = 300;
Formats(end,1).callback = {@set_searchstring,'@ALL','@ALL'};

Prompt{end+1,1} = '';
Formats(end+1,1).type = 'text';

Prompt{end+1,1} = 'Include strings';
Formats(end+1,1).type = 'text';
Prompt(end+1,:) = {'','includestring'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'vector';
Formats(end,1).size = 300;

Prompt{end+1,1} = '';
Formats(end+1,1).type = 'text';

Prompt{end+1,1} = 'Exclude strings';
Formats(end+1,1).type = 'text';
Prompt(end+1,:) = {'','excludestring'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'vector';
Formats(end,1).size = 300;

if doextra == 1
    Prompt{end+1,1} = '';
    Formats(end+1,1).type = 'text';
    
    Prompt(end+1,:) = {'Disregard already processed files','excludeprocessed'};
    Formats(end+1,1).type = 'check';
    
    Prompt(end+1,:) = {'Include ICA results','doICAresult'};
    Formats(end+1,1).type = 'check';
end

[cfg.SEARCH,Cancelled] = inputsdlg(Prompt,'Search settings',Formats,cfg.SEARCH);
if Cancelled == 1
    cfg.SEARCH = [];
    skipprocessing = 1;
else
    pause(0.2);
end

end

function settings = set_searchmode(settings)
    switch settings.searchmode
        case 'ICA'
            settings.searchstring = {'ICA.mat'};
        case 'Filtered'
            settings.searchstring = {'filt.sef'};
        case 'Results'
            settings.searchstring = {'result.'};
        case 'Matrix (.txt)'
            settings.searchstring = {'matrix.txt'};
        case 'Matrix (.mat)'
            settings.searchstring = {'Conn_*.mat'};
        case 'Cartool (.sef)'
            settings.searchstring = {'.sef'};
        case 'Cartool (.ris)'
            settings.searchstring = {'.ris'};
        case 'Cartool (.eph)'
            settings.searchstring = {'.eph'};
        case 'EDF (.edf)'
            settings.searchstring = {'.edf'};
        case 'EGI (.mff)'
            settings.searchstring = {'.mff'};
        case 'EGI (.raw)'
            settings.searchstring = {'.raw'};
        case 'EGI (.mat)'
            settings.searchstring = {'.mat'};
        case 'Text (.txt)'
            settings.searchstring = {'.txt'};
        case 'Elekta (.fif)'
            settings.searchstring = {'.fif';'.fiff'};
        case 'Info-file'
            settings.searchstring = {'.info'};
        case 'Verbose-file'
            settings.searchstring = {'bad.vrb'};
            settings.excludestring = {'ICA_bad.vrb'};
        case 'ICA-Container'
            settings.searchstring = {'ICA.mat'};
        case 'MRI (.hdr)'
            settings.searchstring = {'.hdr'};
    end
    if ~isfield(settings,'searchstring')
        settings.searchstring = {''};
    end
    if ~isfield(settings,'excludestring')
        settings.excludestring = {''};
    end
    if ~isfield(settings,'includestring')
        settings.includestring = {''};
    end
end

function settings = set_searchstring(settings)
    if length(settings.searchstring) == 1
        switch settings.searchstring{1}
            case 'ICA.mat'
                settings.searchmode = 'ICA';
            case 'filt.sef'
                settings.searchmode = 'Filtered';
            case 'result.'
                settings.searchmode = 'Results';
            case 'matrix.txt'
                settings.searchmode = 'Matrix (.txt)';
            case 'Conn_*.mat'
                settings.searchmode = 'Matrix (.mat)';
            case '.sef'
                settings.searchmode = 'Cartool (.sef)';
            case '.ris'
                settings.searchmode = 'Cartool (.ris)';
            case '.edf'
                settings.searchmode = 'EDF (.edf)';
            case '.mff'
                settings.searchmode = 'EGI (.mff)';
            case '.raw'
                settings.searchmode = 'EGI (.raw)';
            case '.mat'
                settings.searchmode = 'EGI (.mat)';
            case '.txt'
                settings.searchmode = 'Text (.txt)';
            case '.fif'
                settings.searchmode = 'Elekta (.fif)';
            case '.fiff'
                settings.searchmode = 'Elekta (.fif)';
            case '_bad.vrb'
                settings.searchmode = 'Verbose-file';
            case '.info'
                settings.searchmode = 'Info-file';
            case '.hdr'
                settings.searchmode = 'MRI (.hdr)';
            otherwise
                settings.searchmode = 'strings';
        end
    else
        settings.searchmode= 'strings';
    end
end