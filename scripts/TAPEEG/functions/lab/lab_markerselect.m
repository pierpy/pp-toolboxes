% Helper script for selecting markers
%
% lab_markerselect(settings,cfg,header)
%
% Written by F. Hatz 2012

function settings = lab_markerselect(settings,cfg,header)

if ~isfield(settings,'markerexclude') | ~isfield(settings,'markerinclude')
    return
end

if isfield(header,'events') & isfield(header.events,'TYP') & ~isempty(header.events.TYP)
    if ~isfield(cfg,'MAIN') | ~isfield(cfg.MAIN,'auto') | cfg.MAIN.auto == 0
        if ~isempty(settings.markerexclude) & strcmp(settings.markerexclude{1,1},'*')
            strlist = unique(header.events.TYP);
            strlist = strlist(:)';
            defaultanswer = [];
            [selection] = listdlg('PromptString','Select Markers:','SelectionMode','multiple', ...
                'ListString',strlist,'InitialValue',defaultanswer);
            if ~isempty(selection)
                settings.markerexclude = strlist(1,selection)';
            else
                settings.markerexclude = [];
            end
            clearvars selection strlist
            pause(0.2);
        end
        if ~isempty(settings.markerinclude) & strcmp(settings.markerinclude{1,1},'*')
            strlist = unique(header.events.TYP);
            strlist = strlist(:)';
            defaultanswer = [];
            [selection] = listdlg('PromptString','Select Markers:','SelectionMode','multiple', ...
                'ListString',strlist,'InitialValue',defaultanswer);
            if ~isempty(selection)
                settings.markerinclude = strlist(1,selection)';
            else
                settings.markerinclude = [];
            end
            clearvars selection strlist
            pause(0.2);
        end
        if ~isempty(settings.markerstart) & strcmp(settings.markerstart,'*')
            strlist = unique(header.events.TYP);
            strlist = strlist(:)';
            defaultanswer = [];
            [selection] = listdlg('PromptString','Select Markers:','SelectionMode','single', ...
                'ListString',strlist,'InitialValue',defaultanswer);
            if ~isempty(selection)
                settings.markerstart = strlist(1,selection)';
            else
                settings.markerstart = [];
            end
            clearvars selection strlist
            pause(0.2);
        end
        if ~isempty(settings.markerstop) & strcmp(settings.markerstop,'*')
            strlist = unique(header.events.TYP);
            strlist = strlist(:)';
            defaultanswer = [];
            [selection] = listdlg('PromptString','Select Markers:','SelectionMode','single', ...
                'ListString',strlist,'InitialValue',defaultanswer);
            if ~isempty(selection)
                settings.markerstop = strlist(1,selection)';
            else
                settings.markerstop = [];
            end
            clearvars selection strlist
            pause(0.2);
        end
    end
    if ~isempty(settings.markerexclude) & strcmp(settings.markerexclude{1,1},'*')
        settings.markerexclude = [];
    end
    if ~isempty(settings.markerinclude) & strcmp(settings.markerinclude{1,1},'*')
        settings.markerinclude = [];
    end
    if ~isempty(settings.markerstop) & strcmp(settings.markerstop,'*')
        settings.markerstop = [];
    end
    if ~isempty(settings.markerstart) & strcmp(settings.markerstart,'*')
        settings.markerstart = [];
    end
    if ~isempty(settings.markerexclude) & strcmp(settings.markerexclude{1,1},'all')
        tmp = unique(header.events.TYP);
        settings.markerexclude = [cellstr('all') tmp(:)'];
    end
    if ~isempty(settings.markerinclude) & strcmp(settings.markerinclude{1,1},'all')
        tmp = unique(header.events.TYP);
        settings.markerinclude = [cellstr('all') tmp(:)'];
    end
end