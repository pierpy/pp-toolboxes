function Microstates = lab_load_microstates(Microstates,header)

global Main_Path

if ~exist('header','var')
    header = [];
end
if isfield(header,'numchannels') & ~isfield(header,'numdatachannels')
    header.numdatachannels = header.numchannels;
end
if ~exist('Microstates','var')
    Microstates = [];
end

List_Micro = {};
if ~isempty(Main_Path) & exist(fullfile(Main_Path,'electrodes'),'dir')
    List_Micro = lab_search(fullfile(Main_Path,'microstates'),{'*.sef','*.eph','*.ris'},true,true,1);
end
List_Micro = cat(1,{'Select File'},List_Micro(:));
List_Micro = cat(1,{''},List_Micro(:));

settings.microstates_file = [];
settings.microstates = Microstates;

Prompt = cell(0,2);
Formats = [];

Prompt(end+1,:) = {'Microstates-File','microstates_file'};
Formats(end+1,1).type = 'list';
Formats(end,1).style = 'popupmenu';
Formats(end,1).format = 'input';
Formats(end,1).items = List_Micro;
Formats(end,1).callback =  {@read_microstates,{'microstates','microstates_file'},'microstates_file',header};
Formats(end,1).size = 300;
Formats(end,1).span = [1 2];

Prompt(end+1,:) = {'Microstates','microstates'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'result';
Formats(end,1).size = [90 20];
Formats(end,1).enable = 'inactive';

Prompt(end+1,:) = {'Reduce',''};
Formats(end+1,1).type = 'button';
Formats(end,1).style = 'pushbutton';
Formats(end,1).size = [80 25];
Formats(end,1).callback = {@reduce_microstates,'microstates','microstates',header};

[settings,Cancelled] = inputsdlg(Prompt,'Load Microstates',Formats,settings);
if Cancelled == 1
    Microstates = [];
    return
end
Microstates = settings.microstates;

    function [microstates,microstates_file] = read_microstates(microstates_file,header)
        if ~isempty(microstates_file) & strcmp(microstates_file,'Select File')
            [microstates_file,microstates_path] = uigetfile({'*.sef;*.eph;*.ris','Microstates-File'},'Select Microstates-File');
            microstates_file = fullfile(microstates_path,microstates_file);
        end
        if isempty(microstates_file) | ~exist(microstates_file,'file')
            microstates_file = '';
            microstates = [];
            return
        end
        microstates = lab_read_data(microstates_file);
        if isempty(microstates)
            return
        end
        if isfield(header,'numdatachannels')
            microstates = reduce_microstates(microstates,header);
        end
    end

    function microstates = reduce_microstates(microstates,header)
        if isfield(header,'numdatachannels') & header.numdatachannels == size(microstates,1)
            return
        end
        channels = cellstr(num2str((1:size(microstates,1))'));
        exclude = lab_get_exclude(size(microstates,1));
        if ~isempty(exclude)
            strdefault = setdiff(1:size(microstates,1),exclude);
        else
            strdefault = 1:size(microstates,1);
        end
        selection = listdlg('promptstring','Included channels:','selectionmode','multiple', ...
            'liststring',channels,'initialvalue',strdefault);
        if ~isempty(selection);
            microstates = microstates(selection,:);
        end
    end
end