function [eeg, gfp_peak_indices, gfp_curve]  = preprocess(self, eeg)
    if (self.confstruct.substract_column_at_start)
        eeg = self.detrendcolumnwise(eeg);
    end   
    gfp_curve = compute_gfp(eeg, self.confstruct.method_GFPeak); %calcola gfp
    [gfp_peak_indices, ~, gfp_curve] = compute_gfp_peaks(gfp_curve, self.confstruct.use_gfp_peaks, self.confstruct.maxima_method, self.confstruct.minDistance, self.confstruct.minAmplitude); %calcolo indici massimo del gfp
    eeg = avgref(eeg, self.confstruct.nch);
    
    if (self.confstruct.setall1)
        eeg = set_gfp_all_1(eeg, gfp_curve);
    end
    if (self.confstruct.normalize)
        eeg = normalize_maps(eeg, 'gfp1', 'GFPL2');
    end
end