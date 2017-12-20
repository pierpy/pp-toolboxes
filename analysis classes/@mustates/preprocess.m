function [eeg, gfp_peak_indices, gfp_curve]  = preprocess(self, eeg, nch)
    if (self.confstruct.substract_column_at_start)
        eeg = self.detrendcolumnwise(eeg);
    end   
    gfp_curve = self.computegfp(eeg, self.confstruct.method_GFPeak); %calcola gfp
    [gfp_peak_indices, ~, gfp_curve] = self.computegfppeaks(gfp_curve); %calcolo indici massimo del gfp
    eeg = self.avgref(eeg, nch);   
    if (self.confstruct.setall1)
        eeg = self.setgfpall1(eeg, gfp_curve);
    end
    if (self.confstruct.normalize)
        eeg = normalizemaps(eeg, 'gfp1', 'GFPL2');
    end
end