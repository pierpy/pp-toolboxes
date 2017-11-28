% [data,header,cfg] = lab_read_fif(filename,cfg) - read Elekta *.fif
% written by F. Hatz 2012

function data = lab_read_info(filename)

data = lab_read_eeginfo(filename);