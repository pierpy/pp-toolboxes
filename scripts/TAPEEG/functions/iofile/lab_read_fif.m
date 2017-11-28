% [data,header,cfg] = lab_read_fif(filename,cfg) - read Elekta *.fif
% written by F. Hatz 2012

function [data,header,cfg] = lab_read_fif(filename,cfg,nodata,segment)

if ~exist('segment','var')
    segment = [];
end
if ~exist('nodata','var')
    nodata = false;
end
if ~exist('cfg','var')
    cfg = [];
end
[data,header,cfg] = lab_read_fiff(filename,cfg,nodata,segment);
