function [m] = M(w,Nchans,ClusterNr)
    m = w*(ClusterNr^(2/Nchans));
end