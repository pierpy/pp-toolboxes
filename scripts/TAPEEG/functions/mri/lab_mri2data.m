function [data,header] = lab_mri2data(header)
    
if isfield(header,'anatomy')
    data = reshape(header.anatomy,[size(header.anatomy,1)*size(header.anatomy,2)*size(header.anatomy,3) size(header.anatomy,4)]);
    if ~isfield(header,'dim')
        header.dim = [size(header.anatomy,1) size(header.anatomy,2) size(header.anatomy,3)];
    end
    header.anatomy = [];
elseif isfield(header,'img')
    data = reshape(header.img,[size(header.img,1)*size(header.img,2)*size(header.img,3) size(header.img,4)]);
    if ~isfield(header,'dim')
        header.dim = [size(header.img,1) size(header.img,2) size(header.img,3)];
    end
    header.img = [];
else
    disp('   Abort: no valid mri')
    data = [];
    return
end

header.channels = char(num2str((1:size(data,1))'));
header.numchannels = size(data,1);
header.numdatachannels = header.numchannels;
header.numauxchannels = 0;
header.samplingrate = [];
header.numtimeframes = size(data,2);
header.EEG_file = 'MRI_file';
header.EEG_filepath = pwd;

Prompt = cell(0,2);
Formats = [];

Prompt(end+1,:) = {[num2str(header.numchannels) ' Channels - ' num2str(header.numtimeframes) ' Timeframes'],''};
Formats(end+1,1).type = 'text';
Formats(end,1).span = [1 2];

Prompt(end+1,:) = {'Filename','EEG_file'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'text';
Formats(end,1).size = 200;

Prompt(end+1,:) = {'Filepath','EEG_filepath'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'dir';
Formats(end,1).size = 200;

Prompt(end+1,:) = {'Samplingrate','samplingrate'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'float';
Formats(end,1).limits = [0 inf];
Formats(end,1).size = 50;

[header,Cancelled] = inputsdlg(Prompt,'Set header',Formats,header);
if Cancelled == 1
    header = [];
else
    pause(0.2);
    [~,~,~,header.EEG_file] = lab_filename(header.EEG_file);
    header.EEG_file = [header.EEG_file '.sef'];
end