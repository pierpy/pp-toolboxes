function eeg = setgfpall1(self, eeg, gfp_curve)
    for i = 1 : size(eeg, 1)
        eeg(i,:) = eeg(i,:)./gfp_curve(i);
    end
    