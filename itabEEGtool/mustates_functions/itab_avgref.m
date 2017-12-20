function eeg = itab_avgref(eeg, nch)
    h = eye(nch)-1/nch;
    eeg = eeg*h;
end