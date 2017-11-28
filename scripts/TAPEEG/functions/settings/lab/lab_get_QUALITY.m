function settings = lab_get_QUALITY(settings,dovector,dochans,doact,doepoch)

if ~exist('doepoch','var')
    doepoch = true;
end
if ~exist('doact','var')
    doact = true;
end
if ~exist('dochans','var')
    dochans = true;
end
if exist('dovector','var') & dovector == 1
    Format = 'vector';
    Size = 100;
else
    Format = 'float';
    Size = 50;
end

if ~exist('settings','var') | isempty(settings)
    settings.Npositiv = 2;
    if dochans == true
        settings.Npositiv = settings.Npositiv + 1;
    end
    if doact == true
        settings.Npositiv = settings.Npositiv + 1;
    end
    settings.cogfreq = [];
    settings.peakfreq = [];
    settings.peakamp = 0.7;
    settings.peak2min = 1.9;
    settings.peakratio = [];
    settings.areapower = 0.8;
    settings.bplast = [];
    if dochans == true
        settings.percentbad = 10;
    else
        settings.percentbad = [];
    end
    settings.badmode = 'interpolated';
    settings.LAPL.lap_maxdistance = 3;
    settings.LAPL.lap_percent = 30;
    if doact == true
        settings.percentbadact = 10;
    else
         settings.percentbadact = [];
    end
    if doepoch == true
        settings.epochquality = [];
    else
        settings.epochquality = [];
    end
end

Prompt = cell(0,3);
Formats = [];

Prompt(end+1,:) = {'Number of positive criteria','Npositiv',''};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'integer';
Formats(end,1).limits = [0 6];
Formats(end,1).size = 50;
Formats(end,1).span = [1 3];

Prompt(end+1,:) = {'  ','',''};
Formats(end+1,1).type = 'text';
Formats(end,1).span = [1 3];

Prompt(end+1,:) = {'Median Frequency >','cogfreq',''};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = Format;
Formats(end,1).limits = [0 9999];
Formats(end,1).size = Size;
Formats(end,1).span = [1 3];

Prompt(end+1,:) = {'Peak Frequency >','peakfreq',''};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = Format;
Formats(end,1).limits = [0 9999];
Formats(end,1).size = Size;
Formats(end,1).span = [1 3];

Prompt(end+1,:) = {'Peak Amplitude >','peakamp',''};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = Format;
Formats(end,1).limits = [0 9999];
Formats(end,1).size = Size;
Formats(end,1).span = [1 3];

Prompt(end+1,:) = {'Peak2Min >','peak2min',''};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = Format;
Formats(end,1).limits = [0 9999];
Formats(end,1).size = Size;
Formats(end,1).span = [1 3];

Prompt(end+1,:) = {'Peak Ratio >','peakratio',''};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = Format;
Formats(end,1).limits = [0 9999];
Formats(end,1).size = Size;
Formats(end,1).span = [1 3];

Prompt(end+1,:) = {'BandPower >','areapower',''};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = Format;
Formats(end,1).limits = [0 9999];
Formats(end,1).size = Size;
Formats(end,1).span = [1 3];

Prompt(end+1,:) = {'BandPower (high) >','bplast',''};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = Format;
Formats(end,1).limits = [0 10];
Formats(end,1).size = Size;
Formats(end,1).span = [1 3];

Prompt(end+1,:) = {'  ','',''};
Formats(end+1,1).type = 'text';
Formats(end,1).span = [1 3];

if dochans == true
    Prompt(end+1,:) = {'Channels','percentbad','%'};
    Formats(end+1,1).type = 'edit';
    Formats(end,1).format = 'float';
    Formats(end,1).limits = [0 100];
    Formats(end,1).size = 50;
    Formats(end,1).callback = {@convert_percent,'percentbad','percentbad'};
    
    Prompt(end+1,:) = {'','badmode',''};
    Formats(end+1,1).type = 'list';
    Formats(end,1).style = 'popupmenu';
    Formats(end,1).format = 'input';
    Formats(end,1).items = {'bad','interpolated'};
    
    Prompt(end+1,:) = {'Laplacian','LAPL',''};
    Formats(end+1,1).type = 'check';
    Formats(end,1).callback = {@set_laplacian,'LAPL','LAPL'};
end

if doact == true
    Prompt(end+1,:) = {'Bad activations','percentbadact','%'};
    Formats(end+1,1).type = 'edit';
    Formats(end,1).format = 'float';
    Formats(end,1).limits = [0 100];
    Formats(end,1).size = 50;
    Formats(end,1).callback = {@convert_percent,'percentbadact','percentbadact'};
    Formats(end,1).span = [1 3];
end

if doepoch == true
    Prompt(end+1,:) = {'Epoch quality','epochquality','%'};
    Formats(end+1,1).type = 'edit';
    Formats(end,1).format = 'float';
    Formats(end,1).limits = [0 100];
    Formats(end,1).size = 50;
    Formats(end,1).callback = {@convert_percent,'epochquality','epochquality'};
    Formats(end,1).span = [1 3];
end

[settings,Cancelled] = inputsdlg(Prompt,'Quality settings',Formats,settings);
if Cancelled == 1
    settings = [];
end

end

function settings = set_laplacian(settings)

Prompt = cell(0,2);
Formats = [];

Prompt(end+1,:) = {'Maximal distance','lap_maxdistance'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'float';
Formats(end,1).limits = [0 100];
Formats(end,1).size = 30;

Prompt(end+1,:) = {'Percent bad/interpolated','lap_percent'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'float';
Formats(end,1).limits = [0 100];
Formats(end,1).size = 30;

[settings,Cancelled] = inputsdlg(Prompt,'Laplacian settings',Formats,settings);
if Cancelled == 1
    settings = [];
end

end

function value = convert_percent(value)
    if value < 1
        value = value * 100;
    end
end