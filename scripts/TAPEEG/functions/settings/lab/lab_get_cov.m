function settings = lab_get_cov(settings,header)

disp ('   Ask for covariance-method')

if ~exist('header','var')
    header = [];
end
if ~exist('settings','var') | ~isfield(settings,'covmethod')
    settings.covmethod = 'Input-File';
end
if ~isfield(settings,'markerexclude')
    settings.markerexclude = {};
end
if ~isfield(settings,'markerinclude')
    settings.markerinclude = {};
end

if isfield(header,'events') & isfield(header.events,'TYP')
    Markers = unique(header.events.TYP);
end

Prompt = cell(0,2);
Formats = [];

Prompt(end+1,:) = {'Covariance matrix','covmethod'};
Formats(end+1,1).type = 'list';
Formats(end,1).style = 'popupmenu';
Formats(end,1).format = 'input';
Formats(end,1).items = {'Input-File','Folder','COV-File'};
Formats(end,1).size = 100;
Formats(end,1).callback = {@set_covmethod,'@ALL','@ALL',header};

Prompt(end+1,:) = {'Range','range'};
Formats(end+1,1).type = 'check';
Formats(end,1).callback = {@get_range,'range','range','covmethod'};

Prompt(end+1,:) = {'File','covfile'};
Formats(end+1,1).type = 'check';
Formats(end,1).callback = {@load_covfile,'covfile','covfile','covmethod',header};

Formats(end+1,1).type = 'none';
Formats(end,1).span = [1 3];

Prompt{end+1,1} = 'Markers to exclude (strings / all)';
Formats(end+1,1).type = 'text';
Formats(end,1).span = [1 3];

Prompt(end+1,:) = {'', 'markerexclude'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'vector';
Formats(end,1).size = 300;
Formats(end,1).limits = [0 1];
if exist('Markers','var')
    Formats(end,1).items = Markers;
end
Formats(end,1).span = [1 3];

Prompt{end+1,1} = 'Markers to include (strings / all)';
Formats(end+1,1).type = 'text';
Formats(end,1).span = [1 3];

Prompt(end+1,:) = {'', 'markerinclude'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'vector';
Formats(end,1).size = 300;
Formats(end,1).limits = [0 1];
if exist('Markers','var')
    Formats(end,1).items = Markers;
end
Formats(end,1).span = [1 3];

[settings,Cancelled] = inputsdlg(Prompt,'Covariance matrix',Formats,settings);
if Cancelled == 1
    settings = [];
end

end

function settings = set_covmethod(settings,header)
   if strcmp(settings.covmethod,'Input-File')
       settings.range = get_range(settings.range,settings.covmethod);
       settings.covfile = [];
   elseif strcmp(settings.covmethod,'COV-File')
       settings.covfile = lab_load_covfile(settings.covfile,header);
       if ~isempty(settings.covfile)
           settings.range = [];
           settings.markerexclude = {};
           settings.markerinclude = {};
       end
   end
end

function range = get_range(range,covmethod)
   if strcmp(covmethod,'Input-File')
       if ~isempty(range)
           tmp.start = range(1);
       else
           tmp.start = 1;
       end
       if length(range) > 1
           tmp.stop = range(2);
       else
           tmp.stop = [];
       end
       
       Prompt = cell(0,2);
       Formats = [];
       
       Prompt(end+1,:) = {'Enter range for covariance-calculation:',''};
       Formats(end+1,1).type = 'text';
       
       Prompt(end+1,:) = {'Start-TF (empty = start):','start'};
       Formats(end+1,1).type = 'edit';
       Formats(end,1).format = 'integer';
       Formats(end,1).limits = [1 inf];
       Formats(end,1).size = 70;
       
       Prompt(end+1,:) = {'Stop-TF (empty = end):  ','stop'};
       Formats(end+1,1).type = 'edit';
       Formats(end,1).format = 'integer';
       Formats(end,1).limits = [1 inf];
       Formats(end,1).size = 70;
       
       [tmp,Cancelled] = inputsdlg(Prompt,'Covariance range',Formats,tmp);
       if Cancelled == 1
           range = [];
           return
       end
       if isempty(tmp.start)
           range = [];
       elseif ~isempty(tmp.stop)
           range = [tmp.start tmp.stop];
       else
           range = tmp.start;
       end
   else
       range = [];
   end
end

function covfile = load_covfile(covfile,covmethod,header)
   if strcmp(covmethod,'COV-File')
       covfile = lab_load_covfile(covfile,header);
   else
       covfile = [];
   end
end