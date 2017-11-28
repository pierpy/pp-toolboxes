function [m] = computem(self, w,Nchans,ClusterNr)
    m = w*(ClusterNr^(2/Nchans));
end