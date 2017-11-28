function [gfp_peak_indices, gfp_peak_values, gfp_curve] = compute_gfp_peaks(gfp_curve, use_gfp_peaks, gfp_peaks_method, minDistance, minAmplitude)
    if use_gfp_peaks
        switch gfp_peaks_method
            case 'simple'
                [gfp_peak_values, gfp_peak_indices] = findpeaks(gfp_curve);

            case 'distance'
               [gfp_peak_values, gfp_peak_indices] = findpeaks(gfp_curve,'MinPeakDistance', minDistance);
            case 'amplitude'
                [gfp_peak_values, gfp_peak_indices] = findpeaks(gfp_curve,'MinPeakHeight', minAmplitude);
            case 'both'
                [gfp_peak_values, gfp_peak_indices] = findpeaks(gfp_curve,'MinPeakDistance', minDistance, 'MinPeakHeight', minAmplitude);
        end        
        
    else
        gfp_peak_indices = 1:length(gfp_curve);
        gfp_peak_values = gfp_curve;
    end
    gfp_curve = gfp_curve;                                