% Script to select multiple eeg-file for sequential evaluation using
% lab_plot_eeg
%
% Written by F. Hatz 2013

function lab_evaluate_eeg

button = questdlg('Search or Select File(s)?','','Search','File','File');

if strcmp(button,'File')
    [Filename,Filepath] = uigetfile('*.*','Select EEG/MEG-file','Multiselect','on');
    if isnumeric(Filename) & Filename == 0
        return
    end
    if iscell(Filename)
        for i = 1:size(Filename,2)
            Filelist{1,i} = fullfile(Filepath,Filename{1,i});
        end
    elseif ischar(Filename)
        Filelist{1,1} = fullfile(Filepath,Filename);
    end
    clearvars Filename Filepath
elseif strcmp(button,'Search')
    tmp = lab_search_files;
    Filelist = tmp.Filelist;
    clearvars tmp
end

cfg = [];
for Nfile = 1:size(Filelist,2)
    [data,header,cfg] = lab_read_data(Filelist{1,Nfile},cfg);
    if ~isempty(data)
        lab_plot_eeg(data,header,cfg,1);
    end
end