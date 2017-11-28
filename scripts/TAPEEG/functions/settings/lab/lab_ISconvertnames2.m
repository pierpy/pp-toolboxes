function settings = lab_ISconvertnames2(settings)
   switch settings.IS_file
        case 'Cartool (.is)'
            settings.IS_file = '.is';
        case 'sLoreta (.spinv)'
            settings.IS_file = '.spinv';
        case 'Leadfield (LF.bin)'
            settings.IS_file = 'LF.bin';
        case 'Headmodel (.mat)'
            settings.IS_file = '.mat';
        case 'MRI-file (.hdr)'
            settings.IS_file = '.hdr';
        case 'MRI-file (.nii)'
            settings.IS_file = '.nii';
        case 'MRI-file (iso.fif)'
            settings.IS_file = 'iso.fif';
        case 'MRI-file (mri.fif)'
            settings.IS_file = 'mri.fif';
        case 'MRI-file (dicom)'
            settings.IS_file = 'dicom';
        case 'Template'
            if exist('IS_file','var') & exist(IS_file,'file')
                settings.IS_file = IS_file;
            end
    end
    if isfield(settings,'SPI_file')
        switch settings.SPI_file
            case 'Solutionpoints (.spi)'
                settings.SPI_file = '.spi';
            case 'MRI-file (IS-file)'
                if settings.IndividualFiles == false
                    settings.SPI_file = lab_convertmriname(settings.IS_file);
                else
                    settings.SPI_file = '.hdr';
                end
            case 'MRI-file (spi.hdr)'
                settings.SPI_file = 'spi.hdr';
            case 'MRI-file (spi.nii)'
                settings.SPI_file = 'spi.nii';
        end
    end
    if isfield(settings,'LOC_file')
        switch settings.LOC_file
            case 'Locations in input file'
                settings.LOC_file = 'data';
            case 'Electrodes-file (.xyz)'
                settings.LOC_file = '.xyz';
            case 'Electrodes-file (.els)'
                settings.LOC_file = '.els';
        end
    end
    if isfield(settings,'MRI_file')
        switch settings.MRI_file
            case 'MRI-file (IS-file)'
                if settings.IndividualFiles == false
                    settings.MRI_file = lab_convertmriname(settings.IS_file);
                else
                    settings.MRI_file = '.hdr';
                end
            case 'MRI-file (.hdr)'
                settings.MRI_file = '.hdr';
        end
    end
    if isfield(settings,'ROIS_file')
        switch settings.ROIS_file
            case 'MRI-Atlas (.nii)'
                settings.ROIS_file = '.nii';
            case 'MRI-Atlas (.hdr)'
                settings.ROIS_file = '.hdr';
            case 'ROIS-file (.rois)'
                settings.ROIS_file = '.rois';
            case 'LOC-file (.xyz)'
                settings.ROIS_file = '.xyz';
        end
    end
end