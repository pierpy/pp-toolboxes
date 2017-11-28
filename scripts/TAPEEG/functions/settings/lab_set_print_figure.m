function settings = lab_set_print_figure(settings)

if ~exist('settings','var') | isempty(settings) | ~isfield(settings,'Filename')
    settings.Filename = '';
    settings.docrop = true;
    settings.resolution = 300;
    settings.antialiasing = 2;
    settings.quality = 70;
    settings.renderer = lower(get(gcf,'Renderer'));
    if ~strcmp(settings.renderer,'opengl')
        settings.renderer = 'painters';
    end
    settings.transparent = true;
    settings.native = false;
end

Prompt = cell(0,2);
Formats = [];

Prompt(end+1,:) = {'Save picture to','Filename'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'file';
Formats(end,1).limits = [1 0];
Formats(end,1).size = [-1 0];
Formats(end,1).span = [1 2];

Formats(end+1,1).type = 'none';
Formats(end,1).span = [1 2];

Prompt(end+1,:) = {'Crop borders','docrop'};
Formats(end+1,1).type = 'check';

Prompt(end+1,:) = {'Quality (value 1-100)','quality'};
Formats(end+1,1).type = 'range';
Formats(end,1).limits = [0 100];

Prompt(end+1,:) = {'Resolution (dpi)','resolution'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'integer';
Formats(end,1).limits = [1 9999];
Formats(end,1).size = 50;

Prompt(end+1,:) = {'Antialiasing','antialiasing'};
Formats(end+1,1).type = 'list';
Formats(end,1).style = 'popupmenu';
Formats(end,1).items = {'1','2','3','4'};

Formats(end+1,1).type = 'none';
Formats(end,1).span = [1 2];

Prompt(end+1,:) = {'Renderer','renderer'};
Formats(end+1,1).type = 'list';
Formats(end,1).style = 'popupmenu';
Formats(end,1).format = 'input';
Formats(end,1).items = {'auto','opengl','zbuffer','painters'};
Formats(end,1).span = [1 2];

Formats(end+1,1).type = 'none';
Formats(end,1).span = [1 2];

Prompt(end+1,:) = {'PDF / EPS / PNG',''};
Formats(end+1,1).type = 'text';
Formats(end,1).span = [1 2];

Prompt(end+1,:) = {'Colorspace','colorspace'};
Formats(end+1,1).type = 'list';
Formats(end,1).style = 'popupmenu';
Formats(end,1).format = 'input';
Formats(end,1).items = {'RGB','CMYK','Gray'};

Prompt(end+1,:) = {'Transparent background','transparent'};
Formats(end+1,1).type = 'check';

Prompt(end+1,:) = {'Add to existing pdf','addpdf'};
Formats(end+1,1).type = 'check';
Formats(end,1).span = [1 2];

Formats(end+1,1).type = 'none';
Formats(end,1).span = [1 2];

Prompt(end+1,:) = {'Images',''};
Formats(end+1,1).type = 'text';
Formats(end,1).span = [1 2];

Prompt(end+1,:) = {'Resolution of displayed image','native'};
Formats(end+1,1).type = 'check';

[settings,Cancelled] = inputsdlg(Prompt,'Print figure',Formats,settings);
if Cancelled == 1
    settings = [];
end

end
