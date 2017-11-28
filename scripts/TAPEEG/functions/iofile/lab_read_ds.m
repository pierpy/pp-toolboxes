% Wrapper to [data,header,cfg] = lab_read_meg4(filename,cfg) - read CTF-MEG
%
% written by F. Hatz 2012

function [data,header,cfg] = lab_read_ds(filename,cfg)

if ~exist('cfg','var')
    cfg = [];
end

[data,header,cfg] = lab_read_meg4(filename,cfg);