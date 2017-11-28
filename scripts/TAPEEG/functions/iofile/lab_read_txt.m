% Read eeg/meg txt-file (ascii)
%
% [data,header,cfg]=lab_read_txt(txt_file,cfg)
% 
% by F. Hatz, Neurology Basel

function [data,header,cfg]=lab_read_txt(txt_file,cfg)

if ~exist('cfg','var')
    cfg = [];
end

if exist('txt_file','var')
    [txt_file,txt_filepath,~,txt_fileS] = lab_filename(txt_file);
    if isempty(txt_filepath)
        txt_filepath=pwd;
    end
else
    [txt_file,txt_filepath]=uigetfile('*.txt','Select txt-file');
    if txt_file == 0
        data = [];
        header = [];
        return
    else
        [~,~,~,txt_fileS] = lab_filename(txt_file);
    end
end

try
    data = importdata(fullfile(txt_filepath,txt_file));
catch %#ok<CTCH>
    data = 'error';
end
if isstruct(data) & isfield(data,'textdata') & size(data.textdata,1) == 1 & ...
    size(data.textdata,2) == 1
     data.textdata = textscan(data.textdata{1,1},'%s');
     data.textdata = data.textdata{1,1};
     if size(data.textdata,1) == size(data.data,2);
         data.textdata = data.textdata';
     end
end
if isstruct(data) & size(data.data,1) > 1 & min(size(data.textdata)) == 1
    if size(data.textdata,2) > size(data.textdata,1)
        header.channels =  char(data.textdata');
        data = data.data';
    else
        header.channels =  char(data.textdata);
        data = data.data;
    end
    if ~exist('cfg','var') | ~isfield(cfg,'TXT') | ~isfield(cfg.TXT,'samplingrate')
        if exist(fullfile(txt_filepath,[txt_fileS '.info']),'file')
            headertmp = lab_read_eeginfo(fullfile(txt_filepath,[txt_fileS '.info']));
            if isfield(headertmp,'samplingrate')
                cfg.TXT.samplingrate = headertmp.samplingrate;
            else
                cfg.TXT.samplingrate = [];
            end
        else
            cfg.TXT.samplingrate = [];
        end
        if ~isfield(cfg,'MAIN') | ~isfield(cfg.MAIN,'auto') | cfg.MAIN.auto == 0
            Prompt = cell(0,2);
            Formats = [];
            Prompt(end+1,:) = {'Samplingrate','samplingrate'};
            Formats(end+1,1).type = 'edit';
            Formats(end,1).format = 'integer';
            Formats(end,1).limits = [0 9999];
            Formats(end,1).size = 40;
            [cfg.TXT,Cancelled] = inputsdlg(Prompt,'Import txt-file',Formats,cfg.TXT);
            if isempty(cfg.TXT) | Cancelled == 1
                cfg.TXT.samplingrate = [];
            else
                pause(0.2);
            end
        end
    end
    header.numtimeframes = size(data,2);
    header.numchannels = size(data,1);
    header.samplingrate = cfg.TXT.samplingrate;
elseif isnumeric(data) & size(data,1) > 1
    if size(data,1) ~= size(data,2)
        data = data';
        if ~exist('cfg','var') | ~isfield(cfg,'TXT') | ~isfield(cfg.TXT,'labels')
            cfg.TXT.headerlines = 0;
            if exist(fullfile(txt_filepath,[txt_fileS '.info']),'file')
                headertmp = lab_read_eeginfo(fullfile(txt_filepath,[txt_fileS '.info']));
                if isfield(headertmp,'samplingrate')
                    cfg.TXT.samplingrate = headertmp.samplingrate;
                else
                    cfg.TXT.samplingrate = [];
                end
            else
                cfg.TXT.samplingrate = [];
            end
            if min(data(2:end,1) - data(1:end-1,1)) > 0
                cfg.TXT.labels = num2str(data(:,1)');
                cfg.TXT.headerlines = 1;
            elseif size(data,1) == 24
                cfg.TXT.labels = {'Fp2 Fp1 F8 F7 F4 F3 A2 A1 T4 T3 C4 C3 T6 T5 P4 P3 O2 O1 Fz Cz Pz T2 T1 ECG'};
            elseif size(data,1) == 23
                cfg.TXT.labels = {'Fp2 Fp1 F8 F7 F4 F3 A2 A1 T4 T3 C4 C3 T6 T5 P4 P3 O2 O1 Fz Cz Pz T2 T1'};
            else
                cfg.TXT.labels = num2str((1:size(data,1)));
            end
            if ~isfield(cfg,'MAIN') | ~isfield(cfg.MAIN,'auto') | cfg.MAIN.auto == 0
                Prompt = cell(0,2);
                Formats = [];
                
                Prompt(end+1,:) = {'Labels','labels'};
                Formats(end+1,1).type = 'edit';
                Formats(end,1).format = 'text';
                Formats(end,1).limits = [0 6];
                Formats(end,1).size = 350;
                
                Prompt(end+1,:) = {'Header lines','headerlines'};
                Formats(end+1,1).type = 'edit';
                Formats(end,1).format = 'integer';
                Formats(end,1).limits = [0 99];
                Formats(end,1).size = 30;
                
                Prompt(end+1,:) = {'Samplingrate','samplingrate'};
                Formats(end+1,1).type = 'edit';
                Formats(end,1).format = 'integer';
                Formats(end,1).limits = [0 9999];
                Formats(end,1).size = 40;
                
                [cfg.TXT,Cancelled] = inputsdlg(Prompt,'Import txt-file',Formats,cfg.TXT);
                if isempty(cfg.TXT) | Cancelled == 1
                    cfg.TXT.samplingrate = [];
                    cfg.TXT.labels = cellstr(num2str((1:size(data,1))'));
                    cfg.TXT.headerlines = 0;
                else
                    pause(0.2);
                end
            end
            tmp = textscan(cfg.TXT.labels,'%s');
            if ~isempty(tmp)
                cfg.TXT.labels = tmp{1,1};
            else
                cfg.TXT.labels = cellstr(num2str((1:size(data,1))'));
            end
            clearvars tmp
        end
        header.channels = char(cfg.TXT.labels);
        header.numchannels = size(data,1);
        header.samplingrate = cfg.TXT.samplingrate;
        if ~isempty(cfg.TXT.headerlines) & cfg.TXT.headerlines > 0
            data = data(:,cfg.TXT.headerlines+1:end);
        end
        header.numtimeframes = size(data,2);
    else
        header.numchannels = size(data,1);
        header.numauxchannels = 0;
        header.datatype = 'matrix';
    end
elseif isstruct(data) & isfield(data,'textdata') & strcmp(data.textdata{1,1}(1:5),'File:') & ...
        strcmp(data.textdata{2,1}(1:6),'Epoch:')
    data = lab_read_bw_matrices(fullfile(txt_filepath,txt_file));
    if ~isempty(data)
        header.datatype = 'matrix';
        header.numchannels = size(data.matrix,1);
        header.numauxchannels = 0;
        header.subjects = data.name;
        data = data.matrix;
    else
        header = [];
    end
elseif ~isempty(cfg)
    disp('    error reading txt file, no valid eeg/meg data')
    data = [];
    header = [];
elseif isstruct(data) & isfield(data,'textdata') & length(data.textdata{1,1}) >= 9 & ...
        strcmp(data.textdata{1,1}(1:9),'Filename:')
    data = lab_import_bw_results(fullfile(txt_filepath,txt_file));
    header = [];
else
    if ischar(data) & strcmp(data,'error')
        try
            data = lab_read_bw_matrices(fullfile(txt_filepath,txt_file));
        catch %#ok<CTCH>
            disp('    error reading txt file')
            data = [];
        end
        header = [];
    else
        datatmp = lab_read_eeginfo(fullfile(txt_filepath,txt_file));
        if ~isempty(datatmp)
            data = datatmp;
        end
        header = [];
    end
end    

return


