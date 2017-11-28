function settings = lab_set_edges_nodes(settings,skipcolor)

% Define colors for edges and nodes
if ~exist('skipcolor','var')
    skipcolor = false;
end
if ~exist('settings','var') | ~isfield(settings,'esize')
    settings.esize = 1;
    settings.ColorModeE = 'color';
    settings.ColorE = [1 0 0];
    settings.radius = 1;
    settings.ColorModeN = 'color';
    settings.ColorN = [0 0 0];
end
Prompt = {};
Formats = {};
Prompt(end+1,:) = {'Edge size (default = 1)', 'esize'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'integer';
Formats(end,1).limits = [0 99];
Formats(end,1).size = 40;

Formats(end+1,1).type = 'none';
Formats(end,1).span = [1 2];

if skipcolor == false
    Prompt(end+1,:) = {'Color','ColorModeE'};
    Formats(end+1,1).type = 'list';
    Formats(end,1).style = 'popupmenu';
    Formats(end,1).format = 'input';
    Formats(end,1).items = {'color','bluered','autumn','bone','colorcube','cool','copper','gray', ...
        'hot','hsv','jet','iautumn','ibone','icolorcube','icool','icopper','igray','ihot','ihsv','ijet'};
    Formats(end,1).callback = {@set_color,'ColorE','ColorModeE','ColorE',[1 0 0]};
    
    Prompt(end+1,:) = {'','ColorE'};
    Formats(end+1,1).type = 'color';
    Formats(end,1).size = 20;
end

Formats(end+1,1).type = 'none';
Formats(end,1).span = [1 5];

Prompt(end+1,:) = {'Node size (default = 1)', 'radius'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'integer';
Formats(end,1).limits = [0 99];
Formats(end,1).size = 40;

if skipcolor == false
    Prompt(end+1,:) = {'Min', 'minval'};
    Formats(end+1,1).type = 'edit';
    Formats(end,1).format = 'float';
    Formats(end,1).limits = [-inf inf];
    Formats(end,1).size = 40;
    
    Prompt(end+1,:) = {'Max', 'maxval'};
    Formats(end+1,1).type = 'edit';
    Formats(end,1).format = 'float';
    Formats(end,1).limits = [-inf inf];
    Formats(end,1).size = 40;
    
    Prompt(end+1,:) = {'Color','ColorModeN'};
    Formats(end+1,1).type = 'list';
    Formats(end,1).style = 'popupmenu';
    Formats(end,1).format = 'input';
    Formats(end,1).items = {'color','bluered','autumn','bone','colorcube','cool','copper', ...
        'gray','hot','hsv','jet','iautumn','ibone','icolorcube','icool','icopper','igray','ihot','ihsv','ijet'};
    Formats(end,1).callback = {@set_color,'ColorN','ColorModeN','ColorN',[0 0 0]};
    
    Prompt(end+1,:) = {'','ColorN'};
    Formats(end+1,1).type = 'color';
    Formats(end,1).size = 20;
end

settings = inputsdlg(Prompt,'Node/Edges',Formats,settings);
settings.radius = settings.radius * 0.016;
if strcmp(settings.ColorModeE,'color')
    settings.ecolor = settings.ColorE;
else
    settings.ecolor = settings.ColorModeE;
end
if skipcolor == false
    if strcmp(settings.ColorModeN,'color')
        settings.mycmap = lab_create_cmap(settings.ColorN);
    else
        settings.mycmap = lab_create_cmap(settings.ColorModeN);
    end
end
settings.add = 0;
settings.plot_file = [];

end

function color = set_color(colormode,color,default)
   if strcmp(colormode,'color')
       if min(color) == 1
           color = default;
       end
   else
       color = [1 1 1];
   end
end