[filename,filepath] = uigetfile('*.*','Get filename');
filename = fullfile(filepath,filename);
clearvars filepath