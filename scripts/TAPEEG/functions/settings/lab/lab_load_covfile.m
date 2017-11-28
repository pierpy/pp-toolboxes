function Covfile = lab_load_covfile(Covfile,header)

if ~exist('header','var')
    header = [];
end
if isfield(header,'numchannels') & ~isfield(header,'numdatachannels')
    header.numdatachannels = header.numchannels;
end
if ~exist('Covfile','var')
    Covfile = [];
end
settings.covfile = Covfile;

Prompt = cell(0,2);
Formats = [];

Prompt(end+1,:) = {'COV-File','covfile_file'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'file';
Formats(end,1).items = {'*.cov','Covariance-File (.cov)'};
Formats(end,1).limits = [0 1]; % single file get
Formats(end,1).callback =  {@read_covfile,'covfile','covfile_file',header};
Formats(end,1).size = 270;
Formats(end,1).span = [1 2];

Prompt(end+1,:) = {'Covariance','covfile'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'result';
Formats(end,1).size = [90 20];
Formats(end,1).enable = 'inactive';
Formats(end,1).callback = {@plot_covfile,[],'covfile'};

Prompt(end+1,:) = {'Reduce',''};
Formats(end+1,1).type = 'button';
Formats(end,1).style = 'pushbutton';
Formats(end,1).size = [100 25];
Formats(end,1).callback = {@reduce_covfile,'covfile','covfile',header};

[settings,Cancelled] = inputsdlg(Prompt,'Load Covfile',Formats,settings);
if Cancelled == 1
    Covfile = [];
    return
end
Covfile = settings.covfile;

    function covfile = read_covfile(covfile_file,header)
        if exist(covfile_file,'file')
            covfile = lab_read_cov(covfile_file);
            if isempty(covfile)
                return
            end
            if isfield(header,'numdatachannels')
                covfile = reduce_covfile(covfile,header);
            end
        else
            covfile = [];
        end
    end

    function covfile = reduce_covfile(covfile,header)
        if isfield(header,'numdatachannels') & header.numdatachannels == size(covfile,1)
            return
        end
        channels = cellstr(num2str((1:size(covfile,1))'));
        exclude = lab_get_exclude(size(covfile,1));
        if ~isempty(exclude)
            strdefault = setdiff(1:size(covfile,1),exclude);
        else
            strdefault = 1:size(covfile,1);
        end
        selection = listdlg('promptstring','Included channels:','selectionmode','multiple', ...
            'liststring',channels,'initialvalue',strdefault);
        if ~isempty(selection);
            covfile = covfile(selection,selection);
        end
    end
    
    function plot_covfile(covfile)
        if ~isempty(covfile)
            lab_plot_matrix(covfile,false,true);
        end
    end
end