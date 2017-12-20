function eeg = avgref(self, eeg, nch)
    h = eye(nch)-1/nch;
    %eeg(gfp_peak_indices,:) = eeg(gfp_peak_indices,:)*h;
    eeg = eeg*h;