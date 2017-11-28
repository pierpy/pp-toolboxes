function [MRI_file] = lab_get_MRI(MRI_file)

if strcmp(MRI_file,'Select File')
    [MRI_file,MRI_filepath] = uigetfile({'*.hdr','MRI-file (.hdr)'},'Select MRI-file');
    MRI_file = fullfile(MRI_filepath,MRI_file);
end