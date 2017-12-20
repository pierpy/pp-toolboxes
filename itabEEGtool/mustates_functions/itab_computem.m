function [m] = itab_computem(w,Nchans,ClusterNr)
    m = w*(ClusterNr^(2/Nchans));
end