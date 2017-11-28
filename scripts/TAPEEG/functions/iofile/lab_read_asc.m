% Read eeg/meg txt-file (ascii)
%
% [data,header,cfg]=lab_read_txt(txt_file,cfg)
% 
% by F. Hatz, Neurology Basel

function [data,header,cfg]=lab_read_asc(txt_file,cfg)

if ~exist('cfg','var')
    cfg = [];
end

if exist('txt_file','var')
    tmp=strfind(txt_file,filesep);
    if ~isempty(tmp)
        txt_filepath=txt_file(1:tmp(end));
        txt_file=txt_file(tmp(end)+1:end);
    else
        txt_filepath=pwd;
    end
else
    [txt_file,txt_filepath]=uigetfile('*.txt','Select txt_data.txt file');
end

data=importdata(fullfile(txt_filepath,txt_file));
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
        disp ('   Ask for samplingrate')
        prompt={'Samplingrate'};
        name='Samplingrate';
        numlines(1,:) = [1 30];
        defaultanswer={''};
        answer=inputdlg(prompt,name,numlines,defaultanswer);
        pause(0.1)
        if size(answer,1) > 0
            cfg.TXT.samplingrate = str2num(answer{1,1});
        else
            cfg.TXT.samplingrate = [];
        end
    end
    header.numtimeframes = size(data,2);
    header.numchannels = size(data,1);
    header.samplingrate = cfg.TXT.samplingrate;
elseif isnumeric(data) & size(data,1) > 1
    data = data';
    if ~exist('cfg','var') | ~isfield(cfg,'TXT') | ~isfield(cfg.TXT,'labels')
        disp ('   Ask for txt settings')
        prompt={'Labels','Samplingrate','Header lines'};
        name='Labels';
        numlines(1,:) = [2 50];
        numlines(2,:) = [1 30];
        numlines(3,:) = [1 30];
        if size(data,1) == 24
            defaultanswer={'Fp2 Fp1 F8  F7  F4  F3  A2  A1  T4  T3  C4  C3  T6  T5  P4  P3  O2  O1  Fz  Cz  Pz  T2  T1  ECG','','0'};
        elseif size(data,1) == 23
            defaultanswer={'Fp2 Fp1 F8  F7  F4  F3  A2  A1  T4  T3  C4  C3  T6  T5  P4  P3  O2  O1  Fz  Cz  Pz  T2  T1','','0'};
        else
            defaultanswer={num2str(1:size(data,1)),'','0'};
        end
        answer=inputdlg(prompt,name,numlines,defaultanswer);
        pause(0.1)
        if size(answer,1) > 0
            tmp = textscan(answer{1,1},'%s');
            if ~isempty(tmp)
                cfg.TXT.labels = tmp{1,1};
            else
                cfg.TXT.labels = cellstr(num2str((1:size(data,1))'));
            end
            cfg.TXT.samplingrate = str2num(answer{2,1}); %#ok<ST2NM>
            cfg.TXT.headerlines = str2num(answer{3,1}); %#ok<ST2NM>
        else
            cfg.TXT.labels = cellstr(num2str((1:size(data,1))'));
            cfg.TXT.samplingrate = [];
            cfg.TXT.headerlines = 0;
        end
    end
    header.channels = char(cfg.TXT.labels);
    header.numchannels = size(data,1);
    header.samplingrate = cfg.TXT.samplingrate;
    if ~isempty(cfg.TXT.headerlines) & cfg.TXT.headerlines > 0
        data = data(:,cfg.TXT.headerlines+1:end);
    end
    header.numtimeframes = size(data,2);
elseif ~isempty(cfg)
    disp('    error reading txt file, no valid eeg/meg data')
    data = [];
    header = [];
elseif isstruct(data) & isfield(data,'textdata') & length(data.textdata{1,1}) >= 9 & ...
        strcmp(data.textdata{1,1}(1:9),'Filename:')
    data = lab_import_bw_results(fullfile(txt_filepath,txt_file));
    header = [];
elseif isstruct(data) & isfield(data,'textdata') & strcmp(data.textdata{1,1}(1:5),'File:') & ...
        strcmp(data.textdata{2,1}(1:6),'Epoch:')
    data = lab_read_bw_matrices(fullfile(txt_filepath,txt_file));
    header = [];
else
    datatmp = lab_read_eeginfo(fullfile(txt_filepath,txt_file));
    if ~isempty(datatmp)
        data = datatmp;
    end
    header = [];
end    

return


