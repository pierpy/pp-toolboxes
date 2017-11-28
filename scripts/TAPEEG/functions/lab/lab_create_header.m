% Create standard header for data matrix (channels x timeframes)
% 
% header = lab_create_header(data)
%
% written by F. Hatz 2012

function header = lab_create_header(data,short)

if ~exist('short','var')
    short = false;
end
if isempty(data)
    header = [];
    return
end
if size(data,1) == 24
    header.channels = char({'Fp2';'Fp1';'F8';'F7';'F4';'F3';'A2';'A1';'T4';'T3';'C4';'C3';'T6';'T5';'P4';'P3';'O2';'O1';'Fz';'Cz';'Pz';'T2';'T1';'ECG'});
elseif size(data,1) == 23
    header.channels = char({'Fp2';'Fp1';'F8';'F7';'F4';'F3';'A2';'A1';'T4';'T3';'C4';'C3';'T6';'T5';'P4';'P3';'O2';'O1';'Fz';'Cz';'Pz';'T2';'T1'});
else
    header.channels = char(num2str((1:size(data,1))'));
end
header.numchannels = size(data,1);
header.numdatachannels = header.numchannels;
header.numauxchannels = 0;
header.samplingrate = [];
header.numtimeframes = size(data,2);
header.EEG_file = 'EEG_file';
header.EEG_filepath = pwd;

if short == true
    header.EEG_filepath = '';
    return
end

disp ('   Ask for header settings')

Prompt = cell(0,2);
Formats = [];

Prompt(end+1,:) = {[num2str(header.numchannels) ' Channels - ' num2str(header.numtimeframes) ' Timeframes'],''};
Formats(end+1,1).type = 'text';
Formats(end,1).span = [1 2];

Prompt(end+1,:) = {'Filename','EEG_file'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'text';
Formats(end,1).size = 200;

if size(data,1) <= 500
    Prompt(end+1,:) = {'Labels','channels'};
    Formats(end+1,1).type = 'edit';
    Formats(end,1).format = 'text';
    Formats(end,1).limits = [0 9];
    Formats(end,1).size = [60 200];
    Formats(end,1).span = [9 1];
else
    Prompt(end+1,:) = {' ',''};
    Formats(end+1,1).type = 'text';
    Formats(end,1).span = [9 1];
end

Prompt(end+1,:) = {'Filepath','EEG_filepath'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'dir';
Formats(end,1).size = 200;

Prompt(end+1,:) = {'Number of auxillary channels','numauxchannels'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'integer';
Formats(end,1).limits = [0 header.numchannels-1];
Formats(end,1).size = 40;

Prompt(end+1,:) = {'Samplingrate','samplingrate'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'float';
Formats(end,1).limits = [0 inf];
Formats(end,1).size = 50;

[header,Cancelled] = inputsdlg(Prompt,'Set header',Formats,header,2);
if Cancelled == 1
    header = [];
else
    pause(0.2);
    [~,~,~,header.EEG_file] = lab_filename(header.EEG_file);
    header.EEG_file = [header.EEG_file '.sef'];
    if size(header.channels,1) ~= header.numchannels;
        header.channels = char(num2str((1:size(data,1))'));
    end
end