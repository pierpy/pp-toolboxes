function [settings,skipprocessing] = lab_set_xls2plot(settings,data,vars,varnames,Filename)

skipprocessing = 0;

if length(vars) < 2
    skipprocessing = 1;
    disp('Abort: need more variables')
    return
end

vars = regexprep(vars,'_',' ');
[~,~,~,Filename] = lab_filename(Filename);
Filename = regexprep(Filename,'_',' ');

if ~exist('settings','var') | ~isfield(settings,'xvar')
    settings.title = Filename;
    settings.mode = 'Boxplot';
    if isfield(settings,'numresults') & settings.numresults == 0
        settings.xvar = 1;
    else
        settings.xvar = length(vars)+1;
    end
    settings.xlog = false;
    settings.yvar = 1;
    settings.yname = vars(1);
    settings.ylog = false;
    settings.statistics = 'None';
    settings.Legend = [];
    settings.Marker = [];
    settings.XequalY = false;
    settings.colormode = 'color';
    settings.fileformat = 'jpg';
else
    settings.xvar = settings.xvar + 1;
end

if ~isfield(settings,'xname')
    settings.xname = xvar2list(settings);
end

Prompt = {};
Formats = {};

Prompt(end+1,:) = {'Title','title'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'text';
Formats(end,1).size = 150;
Formats(end,1).span = [1 4];

Prompt(end+1,:) = {'Mode','mode'};
Formats(end+1,1).type = 'list';
Formats(end,1).style = 'popupmenu';
Formats(end,1).format = 'input';
Formats(end,1).items = {'Boxplot';'Scatter';'Q-Q-Plot';'Q-Q-Plot-Normal'};
Formats(end,1).callback = {@set_mode,'@ALL','@ALL'};

Prompt(end+1,:) = {'Marker','Marker'};
Formats(end+1,1).type = 'check';
Formats(end,1).callback = {@get_Marker,'Marker','Marker','mode'};

Prompt(end+1,:) = {'Legend','Legend'};
Formats(end+1,1).type = 'check';
Formats(end,1).callback = {@get_Legend,'Legend','Legend','mode'};
Formats(end,1).span = [1 2];

Prompt(end+1,:) = {'',[]};
Formats(end+1,1).type = 'text';

Prompt(end+1,:) = {'Xaxis = Yaxis','XequalY'};
Formats(end+1,1).type = 'check';
Formats(end,1).format = 'integer';
Formats(end,1).callback = {@get_XequalY,{'XequalY','Limits'},'XequalY','mode','Limits'};

Prompt(end+1,:) = {'Limits','Limits'};
Formats(end+1,1).type = 'check';
Formats(end,1).callback = {@get_Limits,'@ALL','@ALL',data};
Formats(end,1).span = [1 2];

Prompt(end+1,:) = {'',[]};
Formats(end+1,1).type = 'text';
Formats(end,1).span = [1 4];

Prompt(end+1,:) = {'X-Axis','xvar'};
Formats(end+1,1).type = 'list';
Formats(end,1).style = 'popupmenu';
Formats(end,1).items = cat(1,cellstr('-VarNames-'),vars(:));
Formats(end,1).callback = {@set_xvar,'@ALL','@ALL'};

Prompt(end+1,:) = {'Name','xname'};
Formats(end+1,1).type = 'table';
Formats(end,1).format = 'table';
Formats(end,1).size = [100 100];
Formats(end,1).span = [1 2];

Prompt(end+1,:) = {'Log','xlog'};
Formats(end+1,1).type = 'check';
Formats(end,1).format = 'integer';
Formats(end,1).callback = {@set_xlog,'@ALL','@ALL'};

Prompt(end+1,:) = {'',[]};
Formats(end+1,1).type = 'text';
Formats(end,1).span = [1 4];

Prompt(end+1,:) = {'Y-Axis','yvar'};
Formats(end+1,1).type = 'list';
Formats(end,1).style = 'listbox';
Formats(end,1).limits = [0 4]; % multi-select
Formats(end,1).items = vars;
Formats(end,1).size = [150 160];
Formats(end,1).callback = {@set_yvar,'@ALL','@ALL'};

Prompt(end+1,:) = {'Name','yname'};
Formats(end+1,1).type = 'table';
Formats(end,1).format = 'table';
Formats(end,1).size = [100 160];
Formats(end,1).span = [1 2];

Prompt(end+1,:) = {'Log','ylog'};
Formats(end+1,1).type = 'table';
Formats(end,1).format = 'table';
% Formats(end,1).items = {{'Log'}};
Formats(end,1).size = [80 160];

Prompt(end+1,:) = {'',[]};
Formats(end+1,1).type = 'text';
Formats(end,1).span = [1 4];

Prompt(end+1,:) = {'Statistics','statistics'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'text';
Formats(end,1).enable = 'inactive';
Formats(end,1).size = 100;
Formats(end,1).callback = {@set_statistics,'statistics','@ALL'};
Formats(end,1).span = [1 4];

Prompt(end+1,:) = {'Fileformat for save','fileformat'};
Formats(end+1,1).type = 'list';
Formats(end,1).style = 'popupmenu';
Formats(end,1).format = 'input';
Formats(end,1).items = {'tif','bmp','emf','pdf','eps','png','jpg'};

[settings,Cancelled] = inputsdlg(Prompt,'XLS Plotting',Formats,settings);
if isempty(settings) | Cancelled == 1
    settings = [];
    skipprocessing = 1;
    return
else
    if settings.xvar > 1
        settings.xvar = settings.xvar - 1;
    else
        settings.xvar = [];
    end
end

    function settings = set_mode(settings)
        settings.xname = xvar2list(settings);
        settings.statistics = 'None';
        if strcmp(settings.mode,'Boxplot')
            settings.xlog = false;
            settings.Legend = [];
            settings.Marker = [];
            settings.XequalY = false;
        elseif strcmp(settings.mode,'Q-Q-Plot-Normal')
            settings.yvar = [];
            settings.yname = {};
            settings.ylog = [];
            settings.Legend = [];
            settings.Marker = [];
            settings.XequalY = false;
        elseif strcmp(settings.mode,'Q-Q-Plot')
            if isempty(settings.yvar)
                settings.yvar = 1;
            else
                settings.yvar = settings.yvar(1);
            end
            settings.yname = vars(settings.yvar);
            settings.ylog = false(length(settings.yvar),1);
            settings.Legend = [];
            settings.Marker = [];
            settings.XequalY = false;
        elseif strcmp(settings.mode,'Scatter')
            settings.Legend.border = true;
            settings.Legend.place = 'NorthEast';
            settings.Marker.size = 36;
            settings.Marker.marker = 'o';
            settings.XequalY = true;
        end
    end

    function settings = set_xvar(settings)
        xvar = settings.xvar - 1;
        if xvar > 0
            if ~isempty(xvar)
                settings.xname = xvar2list(settings);
                settings.xlog = false;
                settings.statistics = 'None';
            end
            if ~isempty(intersect(settings.yvar,xvar))
                [~,tmp] = intersect(settings.yvar,xvar);
                tmp = setdiff(1:length(settings.yvar),tmp);
                if ~isempty(tmp)
                    settings.yvar = settings.yvar(tmp);
                    settings.yname = settings.yname(tmp);
                    settings.ylog = settings.ylog(tmp);
                elseif xvar == 1
                    settings.yvar = 2;
                    settings.yname = vars(2);
                    settings.ylog = false;
                else
                    settings.yvar = 1;
                    settings.yname = vars(1);
                    settings.ylog = false;
                end
            end
        else
            settings.xname = {};
        end
    end

    function settings = set_xlog(settings)
        if ~strcmp(settings.mode,'Boxplot')
            if settings.xlog == true
                settings.xlog = false;
            else
                settings.xlog = true;
            end
        else
            settings.xlog = false;
        end
    end

    function settings = set_yvar(settings)
        xvar = settings.xvar - 1;
        if ~isempty(settings.yvar)
            settings.yname = vars(settings.yvar);
            settings.ylog = false(length(settings.yvar),1);
        end
        if xvar > 0
            if ~isempty(intersect(settings.yvar,xvar))
                [~,tmp] = intersect(settings.yvar,settings.xvar);
                tmp = setdiff(1:length(settings.yvar),tmp);
                if ~isempty(tmp)
                    settings.yvar = settings.yvar(tmp);
                    settings.yname = settings.yname(tmp);
                    settings.ylog = settings.ylog(tmp);
                elseif settings.xvar == 1
                    settings.yvar = 2;
                    settings.yname = vars(2);
                    settings.ylog = false;
                else
                    settings.yvar = 1;
                    settings.yname = vars(1);
                    settings.ylog = false;
                end
            end
        end
        if ~isempty(xvar) & xvar > 0
            xname = vars{xvar};
            tmpx = strfind(xname,'_');
            for j = 1:length(settings.yvar)
                yname = vars{settings.yvar(j)};
                tmpy = strfind(yname,'_');
                if ~isempty(tmpx) & ~isempty(tmpy) & strcmp(yname(1:tmpx(end)),xname(1:tmpy(end)))
                    settings.yname{j} = yname(tmpy(end)+1:end);
                end
            end
        end
    end

    function statistics = set_statistics(settings)
        xvar = settings.xvar - 1;
        if xvar > 0
            if strcmp(settings.mode,'Boxplot')
                groupvalid = 1:size(data,2);
                groupvalid = setdiff(groupvalid,find(isnan(data(xvar,:))));
                groupvalid = setdiff(groupvalid,find(isempty(data(xvar,:))));
                varvalid = 1:size(data,2);
                for j = 1:length(settings.yvar)
                    varvalid = setdiff(varvalid,find(isnan(data(settings.yvar(j),:))));
                    varvalid = setdiff(varvalid,find(isempty(data(settings.yvar(j),:))));
                end
                varvalid = intersect(varvalid,groupvalid);
                tmp = unique(data(xvar,varvalid));
                tmp = sort(tmp);
                if length(tmp) == 2
                    items = {'None','T-test','Mann-Whitney-U'};
                else
                    items = {'None','Anova','KruskalWallis'};
                end
            elseif strcmp(settings.mode,'Scatter')
                items = {'None','Pearson','Kendall','Spearman'};
            else
                statistics = 'None';
                return
            end
            statistics = listdlg('ListString',items,'Name','Statistics','SelectionMode','single','ListSize',[150 80]);
            if ~isempty(statistics)
                statistics = items{statistics};
            else
                statistics = 'None';
            end
        else
            statistics = 'None';
            return
        end
    end

    function list = xvar2list(settings)
        xvar = settings.xvar - 1;
        if xvar > 0
            if strcmp(settings.mode,'Boxplot')
                if ~isempty(varnames) & length(varnames) >= xvar & ~isempty(varnames{xvar})
                    list = varnames{xvar}(:);
                else
                    groupvalid = 1:size(data,2);
                    groupvalid = setdiff(groupvalid,find(isnan(data(xvar,:))));
                    groupvalid = setdiff(groupvalid,find(isempty(data(xvar,:))));
                    varvalid = 1:size(data,2);
                    for j = 1:length(settings.yvar)
                        varvalid = setdiff(varvalid,find(isnan(data(settings.yvar(j),:))));
                        varvalid = setdiff(varvalid,find(isempty(data(settings.yvar(j),:))));
                    end
                    varvalid = intersect(varvalid,groupvalid);
                    tmp = unique(data(xvar,varvalid));
                    tmp = sort(tmp);
                    for j = 1:length(tmp);
                        if isnumeric(tmp(j))
                            tmp2 = num2str(tmp(j));
                        elseif iscell(tmp(j))
                            tmp2 = tmp{j};
                        else
                            tmp2 = tmp(j);
                        end
                        list{j,1} = regexprep(tmp2,'_',' '); %#ok<AGROW>
                    end
                end
            else
                xname = vars{xvar};
                if ~isempty(settings.yvar)
                    yname = vars{settings.yvar(1)};
                    tmpx = strfind(xname,'_');
                    tmpy = strfind(yname,'_');
                    if ~isempty(tmpx) & ~isempty(tmpy) & strcmp(yname(1:tmpx(end)),xname(1:tmpy(end)))
                        xname = xname(tmpx(end)+1:end);
                    end
                end
                list{1} = xname;
            end
        else
            list = {};
        end
    end
end

function Marker = get_Marker(Marker,Mode)

if strcmp(Mode,'Boxplot')
    Marker = [];
    return
end

if isempty(Marker)
    Marker.size = 36;
    Marker.marker = 'o';
end

Prompt = cell(0,2);
Formats = [];

Prompt(end+1,:) = {'Marker','marker'};
Formats(end+1,1).type = 'list';
Formats(end,1).style = 'popupmenu';
Formats(end,1).format = 'input';
Formats(end,1).items = {'o';'.';'*';'x'};

Prompt(end+1,:) = {'',[]};
Formats(end+1,1).type = 'text';

Prompt(end+1,:) = {'Size','size'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'float';
Formats(end,1).size = 50;
Formats(end,1).limits = [0 999];

[Marker,Cancelled] = inputsdlg(Prompt,'Marker',Formats,Marker);
if Cancelled == 1
    Marker = [];
end

end

function Legend = get_Legend(Legend,Mode)

if ~strcmp(Mode,'Scatter')
    Legend = [];
    return
end

if isempty(Legend)
    Legend.border = true;
    Legend.place = 'NorthEast';
end

Prompt = cell(0,2);
Formats = [];

Prompt(end+1,:) = {'Place','place'};
Formats(end+1,1).type = 'list';
Formats(end,1).style = 'popupmenu';
Formats(end,1).format = 'input';
Formats(end,1).items = {'NorthWest','NorthEast','SouthWest','SouthEast','NorthEastOutside','EastOutside','SouthEastOutside','Best'};

Prompt(end+1,:) = {'',[]};
Formats(end+1,1).type = 'text';

Prompt(end+1,:) = {'Border','border'};
Formats(end+1,1).type = 'check';

[Legend,Cancelled] = inputsdlg(Prompt,'Legend',Formats,Legend);
if Cancelled == 1
    Legend = [];
end

end

function [XequalY,Limits] = get_XequalY(XequalY,Mode,Limits)

if ~strcmp(Mode,'Scatter')
    XequalY = false;
    return
end

if XequalY == false
    XequalY = true;
else
    XequalY = false;
end

if XequalY == true
    Limits.x1 = min(Limits.x1,Limits.y1);
    Limits.x2 = max(Limits.x2,Limits.y2);
    Limits.y1 = Limits.x1;
    Limits.y2 = Limits.x2;
end

end

function settings = get_Limits(settings,data)

xvar = settings.xvar - 1;
if xvar > 0
    if strcmp(settings.mode,'Scatter') | ~isfield(settings.Limits,'x1') | isempty(settings.Limits.x1)
        tmp = data(xvar,:);
        settings.Limits.x1 = round(min(tmp(:))*100) / 100;
        settings.Limits.x2 = round(max(tmp(:))*100) / 100;
    end
else
    settings.Limits.x1 = [];
    settings.Limits.x2 = [];
end
if strcmp(settings.mode,'Boxplot') | strcmp(settings.mode,'Scatter')
    tmp = data(settings.yvar,:);
    settings.Limits.y1 = round(min(tmp(:))*100) / 100;
    settings.Limits.y2 = round(max(tmp(:))*100) / 100;
end
if settings.XequalY == true & strcmp(settings.mode,'Scatter')
    settings.Limits.x1 = min(settings.Limits.x1,settings.Limits.y1);
    settings.Limits.x2 = max(settings.Limits.x2,settings.Limits.y2);
    settings.Limits.y1 = settings.Limits.x1;
    settings.Limits.y2 = settings.Limits.x2;
end

Prompt = cell(0,2);
Formats = [];
if strcmp(settings.mode,'Scatter')
    Prompt(end+1,:) = {'X-Limit Low', 'x1'};
    Formats(end+1,1).type = 'edit';
    Formats(end,1).format = 'float';
    Formats(end,1).limits = [-inf inf];
    Formats(end,1).size = 30;
    
    Prompt(end+1,:) = {'X-Limit High', 'x2'};
    Formats(end+1,1).type = 'edit';
    Formats(end,1).format = 'float';
    Formats(end,1).limits = [-inf inf];
    Formats(end,1).size = 30;
end
if strcmp(settings.mode,'Boxplot') | strcmp(settings.mode,'Scatter')
    Prompt(end+1,:) = {'Y-Limit Low', 'y1'};
    Formats(end+1,1).type = 'edit';
    Formats(end,1).format = 'float';
    Formats(end,1).limits = [-inf inf];
    Formats(end,1).size = 30;
    
    Prompt(end+1,:) = {'Y-Limit High', 'y2'};
    Formats(end+1,1).type = 'edit';
    Formats(end,1).format = 'float';
    Formats(end,1).limits = [-inf inf];
    Formats(end,1).size = 30;
end
if ~isempty(Formats)
    [settings.Limits,Cancelled] = inputsdlg(Prompt,'Limits',Formats,settings.Limits,2);
    if Cancelled == 1
        settings.Limits = [];
    end
end

end