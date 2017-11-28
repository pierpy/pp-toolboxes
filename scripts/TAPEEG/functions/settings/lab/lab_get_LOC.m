function [LOC_file] = lab_get_LOC(LOC_file)

if strcmp(LOC_file,'Select File')
    [LOC_file,LOC_filepath] = uigetfile({'*.xyz;*.els','Electrodes-file (.els/.xyz)'}, ...
        'Select LOC-file');
    LOC_file = fullfile(LOC_filepath,LOC_file);
end