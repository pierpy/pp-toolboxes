function [eeg, gfp_peak_indices, gfp_curve]  = itab_preprocess(eeg, nch)
    if (confstruct.substract_column_at_start)
        eeg = itab_detrendcolumnwise(eeg);
    end   
    gfp_curve = itab_computegfp(eeg, confstruct.method_GFPeak); %calcola gfp
    [gfp_peak_indices, ~, gfp_curve] = itab_computegfppeaks(gfp_curve); %calcolo indici massimo del gfp
    eeg = itab_avgref(eeg, nch);   
    if (confstruct.setall1)
        eeg = itab_setgfpall1(eeg, gfp_curve);
    end
    if (confstruct.normalize)
        eeg = itab_normalizemaps(eeg, 'gfp1', 'GFPL2');
    end
end