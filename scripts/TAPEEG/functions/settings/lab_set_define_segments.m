function [cfg,skipprocessing] = lab_set_define_segments(cfg,header)

skipprocessing = 0;

global THeader
if ~isempty(THeader)
    header = THeader;
elseif ~exist('header','var')
    header= [];
end
if ~exist('cfg','var')
    cfg = [];
end
if ~isfield(cfg,'SEG')
    cfg.SEG = [];
end

Prompt = cell(0,2);
Formats = [];

Prompt(end+1,:) = {'Method','select'};
Formats(end+1,1).type = 'list';
Formats(end,1).style = 'listbox';
Formats(end,1).format = 'input';
Formats(end,1).size = [170 90];
Formats(end,1).items = {'Select complete file';'Select by Markers';'Select by Start-Stop-Marker';'Select by EDF-file';'Automated'};
Formats(end,1).callback = {@set_selection,'@ALL','@ALL',header,cfg};
Formats(end,1).span = [4 1];

Prompt(end+1,:) = {'1. Marker','markerstart'};
Formats(end+1,1).type = 'check';
Formats(end,1).callback = {@set_marker,'@ALL','@ALL',header};

Prompt(end+1,:) = {'2. Marker','markerstop'};
Formats(end+1,1).type = 'check';
Formats(end,1).callback = {@set_marker,'@ALL','@ALL',header};

Prompt(end+1,:) = {'',''};
Formats(end+1,1).type = 'text';

Prompt(end+1,:) = {'EDF-File','EDF'};
Formats(end+1,1).type = 'check';
Formats(end,1).callback = {@set_edf,'@ALL','@ALL',header,cfg};
Formats(end,1).span = [1 3];

Prompt(end+1,:) = {'Automated','AUTO'};
Formats(end+1,1).type = 'check';
Formats(end,1).callback = {@lab_set_segments_auto,'AUTO','AUTO',header,cfg};

[cfg.SEG,Cancelled] = inputsdlg(Prompt,'Segments setting',Formats,cfg.SEG,4);
pause(0.2);
if Cancelled == 1
    skipprocessing = 1;
    cfg.SEG = [];
    cfg.SEG.select = 'Select complete file';
    return
end

end

function settings = set_selection(settings,header,cfg)
    if strcmp(settings.select,'Select by Markers') | strcmp(settings.select,'Select by Start-Stop-Marker')
        settings = set_marker(settings,header);
    elseif strcmp(settings.select,'Select by EDF-file')
        settings = set_edf(settings,header,cfg);
    elseif strcmp(settings.select,'Automated')
        settings = set_auto(settings,header,cfg);
    end
end

function settings = set_marker(settings,header)
   if strcmp(settings.select,'Select by Markers') | strcmp(settings.select,'Select by Start-Stop-Marker')
       settings.segments = [];
       settings.EDF = [];
       settings.AUTO = [];
       
       Prompt = cell(0,2);
       Formats = [];
       
       if isfield(settings,'MARK') & isfield(settings.MARK,'edit') & ~isempty(settings.MARK.edit) & ...
               size(settings.MARK.edit,2) >= 7
           markerlist = unique(settings.MARK.edit(:,7));
       else
           markerlist = {};
       end
       if isfield(header,'events') & isfield(header.events,'TYP') & ~isempty(header.events.TYP)
           markerlist = union(markerlist,unique(header.events.TYP));
       end
       if ~isempty(markerlist)
           if strcmp(settings.select,'Select by Start-Stop-Marker')
               markerlist = [cellstr('**edit**') cellstr('*start*') markerlist(:)' cellstr('*end*')];
           else
               markerlist = [cellstr('') markerlist(:)'];
           end
       else
           clearvars markerlist
       end
       if strcmp(settings.select,'Select by Start-Stop-Marker')
           Prompt(end+1,:) = {'Start Marker','markerstart'};
       else
           Prompt(end+1,:) = {'Markers to be included in segments (duration of markers must be > 1)',''};
           Formats(end+1,1).type = 'text';
           Formats(end,1).span = [1 2];
           
           Prompt(end+1,:) = {'1. Marker','markerstart'};
       end
       if exist('markerlist','var')
           Formats(end+1,1).type = 'list';
           Formats(end,1).style = 'popupmenu';
           Formats(end,1).format = 'input';
           Formats(end,1).items = markerlist;
           Formats(end,1).callback = {@lab_get_marker,'markerstart','markerstart'};
       else
           Formats(end+1,1).type = 'edit';
           Formats(end,1).format = 'text';
           Formats(end,1).size = 100;
       end
       
       if strcmp(settings.select,'Select by Start-Stop-Marker')
           Prompt(end+1,:) = {'Stop Marker','markerstop'};
       else
           Prompt(end+1,:) = {'2. Marker','markerstop'};
       end
       if exist('markerlist','var')
           Formats(end+1,1).type = 'list';
           Formats(end,1).style = 'popupmenu';
           Formats(end,1).format = 'input';
           Formats(end,1).items = markerlist;
           Formats(end,1).callback = {@lab_get_marker,'markerstop','markerstop'};
       else
           Formats(end+1,1).type = 'edit';
           Formats(end,1).format = 'text';
           Formats(end,1).size = 100;
       end
       
       [settings,Cancelled] = inputsdlg(Prompt,'Marker setting',Formats,settings,2);
       if Cancelled == 1
           settings.markerstart = [];
           settings.markerstop = [];
       end
   elseif ~strcmp(settings.select,'Select by EDF-file')
       settings.markerstart = [];
       settings.markerstop = [];
   end
end

function settings = set_edf(settings,header,cfg)
   if strcmp(settings.select,'Select by EDF-file')
       settings.segments = [];
       settings.AUTO = [];
       
       if isfield(header,'events') & isfield(header.events,'edfimport')
           domarker = 1;
       else
           domarker = 2;
       end
       if ~isfield(settings,'eegsource') | isempty(settings.eegsource)
           settings.eegsource = ' ';
       end
       settings.EDF.markerstart = settings.markerstart;
       settings.EDF.markerstop = settings.markerstop;
       
       if isfield(settings,'MARK') & isfield(settings.MARK,'edit') & ~isempty(settings.MARK.edit) & ...
               size(settings.MARK.edit,2) >= 7
           markerlist = unique(settings.MARK.edit(:,7));
       else
           markerlist = {};
       end
       if isfield(header,'events') & isfield(header.events,'TYP') & ~isempty(header.events.TYP)
           markerlist = union(markerlist,unique(header.events.TYP));
       end
       if ~isempty(markerlist)
           markerlist = [cellstr('**edit**') cellstr('*start*') markerlist(:)' cellstr('*end*')];
       else
           clearvars markerlist
       end
       
       Prompt = cell(0,2);
       Formats = [];
       
       Prompt(end+1,:) = {'Start Marker','markerstart'};
       if domarker == 1 & exist('markerlist','var')
           Formats(end+1,1).type = 'list';
           Formats(end,1).style = 'popupmenu';
           Formats(end,1).format = 'input';
           Formats(end,1).items = markerlist;
           Formats(end,1).callback = {@lab_get_marker,'markerstart','markerstart'};
       else
           Formats(end+1,1).type = 'edit';
           Formats(end,1).format = 'text';
           Formats(end,1).size = 100;
       end
       
       Prompt(end+1,:) = {'Stop Marker','markerstop'};
       if domarker == 1 & exist('markerlist','var')
           Formats(end+1,1).type = 'list';
           Formats(end,1).style = 'popupmenu';
           Formats(end,1).format = 'input';
           Formats(end,1).items = markerlist;
           Formats(end,1).callback = {@lab_get_marker,'markerstop','markerstop'};
       else
           Formats(end+1,1).type = 'edit';
           Formats(end,1).format = 'text';
           Formats(end,1).size = 100;
       end
        
       Prompt(end+1,:) = {' ',''};
       Formats(end+1,1).type = 'text';
       Formats(end,1).span = [1 3];
       
       Prompt(end+1,:) = {'Settings for creating EDF (if not available):',''};
       Formats(end+1,1).type = 'text';
       Formats(end,1).span = [1 3];
       
       Prompt(end+1,:) = {'EDF-Reference','eegsource'};
       Formats(end+1,1).type = 'list';
       Formats(end,1).style = 'popupmenu';
       Formats(end,1).format = 'input';
       if ~max(strcmp({'channels','mean','median','laplacian','montage','input'},settings.eegsource))
           Formats(end,1).items = {settings.eegsource,'channels','mean','median','laplacian','montage','input'};
       else
           Formats(end,1).items = {'channels','mean','median','laplacian','montage','input'};
       end
       Formats(end,1).callback = {@lab_get_eegsource,{'@ALL'},'@ALL',cfg,header,'MontageEDF.xls'};
       
       Prompt(end+1,:) = {'Montage','montage'};
       Formats(end+1,1).type = 'check';
       Formats(end,1).callback = {@lab_load_montage,'montage','montage',cfg,header,'MontageEDF.xls'};
       
       Prompt(end+1,:) = {'Laplacian','LAPL'};
       Formats(end+1,1).type = 'check';
       Formats(end,1).callback = {@lab_get_laplacian,'LAPL','LAPL'};
       
       [settings.EDF,Cancelled] = inputsdlg(Prompt,'EDF setting',Formats,settings.EDF,3);
       if Cancelled == 1
           settings.eegsource = [];
           settings.montage = [];
           settings.REF = [];
           settings.markerstart = [];
           settings.markerstop = [];
       else
           settings.markerstart = settings.EDF.markerstart;
           settings.markerstop = settings.EDF.markerstop;
       end
   else
       settings.EDF = [];
   end
end

function settings = set_auto(settings,header,cfg)
    if strcmp(settings.select,'Automated')
        settings.EDF = [];
        settings.segments = [];
        settings.markerstop = [];
        settings.markerstart = [];
        if ~isfield(settings,'AUTO')
            settings.AUTO = [];
        end
        settings.AUTO = lab_set_segments_auto(settings.AUTO,header,cfg);
        return
    else
        settings.AUTO = [];
    end
end