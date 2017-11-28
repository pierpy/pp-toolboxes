% written by F. Hatz 2013

function Exclude = lab_load_exclude(Exclude,cfg,header,nodiag)

if isfield(cfg,'MAIN') & isfield(cfg.MAIN,'auto') & cfg.MAIN.auto == 1
    nodiag = 1;
elseif ~exist('nodiag','var')
    nodiag = 0;
end

if ~exist('header','var')
    header = [];
end
if ~exist('cfg','var')
    cfg = [];
end
if ~exist('Exclude','var')
    Exclude = [];
end
if isfield(header,'numdatachannels')
    numchans = header.numdatachannels;
elseif isfield(header,'numchannels')
    numchans = header.numchannels;
elseif isfield(cfg,'EXTRA') & isfield(cfg.EXTRA,'numdatachans')
    numchans = cfg.EXTRA.numdatachans;
elseif isfield(header,'ref_chan')
    numchans = header.ref_chan;
elseif isfield(header,'locs') & isfield(header.locs,'x')
    numchans = length(header.locs.x);
else
    numchans = [];
end
if exist('header','var') & isfield(header,'channels')
    list = cellstr(header.channels);
end

if isempty(Exclude)
    if isfield(cfg,'EXTRA') & isfield(cfg.EXTRA,'exclude') & ~isempty(cfg.EXTRA.exclude) & ...
            (~isfield(header,'numchannels') | max(cfg.EXTRA.exclude) <= header.numchannels)
        Exclude = cfg.EXTRA.exclude;
    elseif ~isempty(numchans)
        Exclude = lab_get_exclude(numchans);
    end
end
if ~isempty(Exclude)
    Exclude = Exclude(:)';
end

if nodiag == 0
    settings.Exclude = Exclude;
    if exist('Exclude_file','var')
        settings.Exclude_file = Exclude_file;
    else
        settings.Exclude_file = [];
    end
    Prompt = cell(0,2);
    Formats = [];
    
    Prompt(end+1,:) = {'Excluded Channels - File','Exclude_file'};
    Formats(end+1,1).type = 'edit';
    Formats(end,1).format = 'file';
    Formats(end,1).items = {'*.xls;*.xlsx','Excluded-File (.xls,.xlsx)'};
    Formats(end,1).limits = [0 1]; % single file get
    Formats(end,1).callback =  {@read_exclude,'Exclude','Exclude_file','Exclude',numchans};
    Formats(end,1).size = 270;
    
    Prompt(end+1,:) = {'Channels','Exclude'};
    Formats(end+1,1).type = 'edit';
    Formats(end,1).format = 'vector';
    Formats(end,1).size = 340;
    Formats(end,1).limits = [-inf inf];
    if exist('list','var')
        Formats(end,1).items = list;
    end
    
    Prompt(end+1,:) = {'Save',''};
    Formats(end+1,1).type = 'button';
    Formats(end,1).style = 'pushbutton';
    Formats(end,1).size = [100 25];
    Formats(end,1).callback = {@save_exclude,'Exclude_file','Exclude',numchans};
    
    [settings,Cancelled] = inputsdlg(Prompt,'Excluded Channels',Formats,settings);
    if Cancelled == 1
        Exclude = [];
        return
    end
    Exclude = settings.Exclude;
end

   function Exclude = read_exclude(Exclude_file,Exclude,numchans)
        if ~isempty(Exclude_file) & exist(Exclude_file,'file')
            if ispc
                [Exclude,~,xlsinput] = xlsread(Exclude_file,1);
            else
                [Exclude,~,xlsinput] = xlsread(Exclude_file,1,'','basic');
            end
            if size(xlsinput,2) == 2 & (~exist('numchans','var') | isempty(numchans))
                Exclude = str2num(xlsinput{1,2}); %#ok<ST2NM>
                if ischar(Exclude)
                    Exclude = str2num(Exclude); %#ok<ST2NM>
                end
            elseif size(xlsinput,2) == 2 & ~isempty(find(Exclude==numchans,1))
                Exclude = xlsinput{find(Exclude==numchans,1),2};
                if ischar(Exclude)
                    Exclude = str2num(Exclude); %#ok<ST2NM>
                end
            elseif size(xlsinput,2) == 2 & ischar(xlsinput{1,2})
                Exclude= [];
            elseif size(Exclude,1) > 1 & size(Exclude,2) == 1
                Exclude = Exclude';
            elseif size(Exclude,1) > 1
                Exclude = Exclude(1,:);
            end
       end
   end
   
   function Exclude_file = save_exclude(Exclude,numchans)
       if isempty(Exclude)
           return
       end
       [Exclude_file,Exclude_filepath] = uiputfile('*.xls;*.xlsx','Select File to store');
       if Exclude_file ~= 0
           Exclude_file = fullfile(Exclude_filepath,Exclude_file);
           if exist('numchans','var') & ~isempty(numchans)
               xlsout{1,1} = numchans;
               xlsout{1,2} = num2str(Exclude);
           else
               xlsout = Exclude;
           end
           lab_write_xls(Exclude_file,xlsout);
       else
           Exclude_file = '';
       end
   end
end