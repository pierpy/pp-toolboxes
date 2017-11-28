function settings = lab_ISconvertnames(settings)
   switch settings.IS_file
        case '.is'
            settings.IS_file = 'Cartool (.is)';
        case '.spinv'
            settings.IS_file = 'sLoreta (.spinv)';
        case 'LF.bin'
            settings.IS_file = 'Leadfield (LF.bin)';
        case '.mat'
            settings.IS_file = 'Headmodel (.mat)';
        case '.hdr'
            settings.IS_file = 'MRI-file (.hdr)';
        case '.nii'
            settings.IS_file = 'MRI-file (.nii)';
        case 'iso.fif'
            settings.IS_file = 'MRI-file (iso.fif)';
        case 'mri.fif'
            settings.IS_file = 'MRI-file (mri.fif)';
        case 'dicom'
            settings.IS_file = 'MRI-file (dicom)';
    end
    if isfield(settings,'SPI_file')
        switch settings.SPI_file
            case '.spi'
                settings.SPI_file = 'Solutionpoints (.spi)';
            case lab_convertmriname(settings.IS_file)
                settings.SPI_file = 'MRI-file (IS-file)';
            case '.hdr'
                settings.SPI_file = 'MRI-file (IS-file)';
            case 'spi.hdr'
                settings.SPI_file = 'MRI-file (spi.hdr)';
            case 'spi.nii'
                settings.SPI_file = 'MRI-file (spi.nii)';
        end
    end
    if isfield(settings,'LOC_file')
        switch settings.LOC_file
            case 'data'
                settings.LOC_file = 'Locations in input file';
            case '.xyz'
                settings.LOC_file = 'Electrodes-file (.xyz)';
            case '.els'
                settings.LOC_file = 'Electrodes-file (.els)';
        end
    end
    if isfield(settings,'MRI_file')
        switch settings.MRI_file
            case lab_convertmriname(settings.IS_file)
                settings.MRI_file = 'MRI-file (IS-file)';
            case '.hdr'
                if strcmp(settings.IS_file,'MRI-file (.hdr)')
                    settings.MRI_file = 'MRI-file (IS-file)';
                else
                    settings.MRI_file = 'MRI-file (.hdr)';
                end
        end
    end
    if isfield(settings,'ROIS_file')
        switch settings.ROIS_file
            case '.nii'
                settings.ROIS_file = 'MRI-Atlas (.nii)';
            case '.hdr'
                settings.ROIS_file = 'MRI-Atlas (.hdr)';
            case '.rois'
                settings.ROIS_file = 'ROIS-file (.rois)';
            case '.xyz'
                settings.ROIS_file = 'LOC-file (.xyz)';
            case '.els'
                settings.ROIS_file = 'LOC-file (.els)';
        end
    end
end