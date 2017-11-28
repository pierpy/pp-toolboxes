function [eeg, gfp_peak_indices, gfp_curve] = modmaps_preprocess(eeg, confstruct)

    if (confstruct.substract_column_at_start)
        eeg = substract_column_mean_epochwise(eeg, confstruct.tf);
    end
    
    gfp_curve = compute_gfp(eeg, confstruct.method_GFPeak); %calcola gfp
    [gfp_peak_indices, ~, gfp_curve] = compute_gfp_peaks(gfp_curve, confstruct.use_gfp_peaks, confstruct.maxima_method, confstruct.minDistance, confstruct.minAmplitude); %calcolo indici massimo del gfp
    eeg = avgref(eeg, confstruct.nch);
    
    if (confstruct.setall1)
        eeg = set_gfp_all_1(eeg, gfp_curve);
    end
    if (confstruct.normalize)
        eeg = normalize_maps(eeg, 'gfp1', 'GFPL2');
    end