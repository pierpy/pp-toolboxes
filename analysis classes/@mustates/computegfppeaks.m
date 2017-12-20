function [gfp_peak_indices, gfp_peak_values, gfp_curve] = computegfppeaks(self, gfp_curve)
    if self.confstruct.use_gfp_peaks
        switch self.confstruct.maxima_method
            case 'simple'
                [gfp_peak_values, gfp_peak_indices] = findpeaks(gfp_curve);
            case 'distance'
                [gfp_peak_values, gfp_peak_indices] = findpeaks(gfp_curve,'MinPeakDistance', self.confstruct.minDistance);
            case 'amplitude'
                [gfp_peak_values, gfp_peak_indices] = findpeaks(gfp_curve,'MinPeakHeight', self.confstruct.minAmplitude);
            case 'both'
                [gfp_peak_values, gfp_peak_indices] = findpeaks(gfp_curve,'MinPeakDistance', self.confstruct.minDistance, 'MinPeakHeight', self.confstruct.minAmplitude);
        end           
    else
        gfp_peak_indices = 1:length(gfp_curve);
        gfp_peak_values = gfp_curve;
    end
    gfp_curve = gfp_curve;                                