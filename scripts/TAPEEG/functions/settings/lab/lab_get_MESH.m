function settings = lab_get_MESH(settings,method)

if ~exist('method','var')
    method = '';
end

if ~exist('settings','var') | ~isfield(settings,'Nscalp')
   settings.Nscalp = 3000;
   settings.Nskull = 3000;
   settings.Nbrain = 3000;
end
if ~exist('settings','var') | ~isfield(settings,'SkullDist')
    settings.SkullDist = 4;
end
if isfield(settings,'tissue') & ~isempty(settings.tissue)
    settings.SkullNr = find(strcmp(settings.tissue,'skull'));
elseif ~isfield(settings,'SkullNr')
    settings.SkullNr = 0;
end
if ~isfield(settings,'OtherDist')
    settings.OtherDist = 1;
end

Prompt = cell(0,2);
Formats = [];

if ~strcmp(method,'BEM') & ~strcmp(method,'FEM-simbio')
    Prompt(end+1,:) = {'Number of vertices',''};
    Formats(end+1,1).type = 'text';
    
    Prompt(end+1,:) = {'Scalp','Nscalp'};
    Formats(end+1,1).type = 'edit';
    Formats(end,1).format = 'integer';
    Formats(end,1).limits = [0 999999];
    Formats(end,1).size = 50;
    
    Prompt(end+1,:) = {'Skull','Nskull'};
    Formats(end+1,1).type = 'edit';
    Formats(end,1).format = 'integer';
    Formats(end,1).limits = [0 999999];
    Formats(end,1).size = 50;
    
    Prompt(end+1,:) = {'Brain','Nbrain'};
    Formats(end+1,1).type = 'edit';
    Formats(end,1).format = 'integer';
    Formats(end,1).limits = [0 999999];
    Formats(end,1).size = 50;
end

if strcmp(method,'openmeeg') | strcmp(method,'BEM')
    Formats(end+1,1).type = 'none';
    
    Prompt(end+1,:) = {'BEM-Mesh (openmeeg)',''};
    Formats(end+1,1).type = 'text';
    
    Prompt(end+1,:) = {'Min thickness skull compartment (mm)','SkullDist'};
    Formats(end+1,1).type = 'edit';
    Formats(end,1).format = 'float';
    Formats(end,1).limits = [0 inf];
    Formats(end,1).size = 50;
    
    Prompt(end+1,:) = {'Number of skull mesh (0 = auto)','SkullNr'};
    Formats(end+1,1).type = 'edit';
    Formats(end,1).format = 'integer';
    Formats(end,1).limits = [0 10];
    Formats(end,1).size = 50;
    
    Prompt(end+1,:) = {'Min thickness other compartments (mm)','OtherDist'};
    Formats(end+1,1).type = 'edit';
    Formats(end,1).format = 'float';
    Formats(end,1).limits = [0 inf];
    Formats(end,1).size = 50;
end

if strcmp(method,'FEM-simbio')
    Prompt(end+1,:) = {'Hexahedral-Mesh (FEM-simbio)',''};
    Formats(end+1,1).type = 'text';
    
    Prompt(end+1,:) = {'Shift','shift'};
    Formats(end+1,1).type = 'edit';
    Formats(end,1).format = 'float';
    Formats(end,1).limits = [0 inf];
    Formats(end,1).size = 50;
    
    Prompt(end+1,:) = {'Resolution','resolution'};
    Formats(end+1,1).type = 'edit';
    Formats(end,1).format = 'integer';
    Formats(end,1).limits = [1 10];
    Formats(end,1).size = 50;
end

[settings,Cancelled] = inputsdlg(Prompt,'Mesh settings',Formats,settings);
if Cancelled == 1
    settings = [];
end