% Preprocess coregistration for  pre/postprocessing in TAPEEG, used by
% lab_tapeeg
%
% written by F. Hatz 2012

function lab_prepare_coreg(Filelist,settings)

for filenr = 1:size(Filelist,2)
    [~,header] = lab_read_data(Filelist{1,filenr});
    
    % set is-files for calculation
    IS_file = settings.IS_file;
    SPI_file = settings.SPI_file;
    ROIS_file = settings.ROIS_file;
    MRI_file = settings.MRI_file;
    
    % look for patient IS-files
    tmp=strfind(settings.IS_file,filesep);
    if isempty(tmp) & isfield(settings,'IndividualFiles') & settings.IndividualFiles == true
        IS_file = lab_find_individualIS(IS_file,SPI_file,ROIS_file,MRI_file,header.EEG_filepath);
    end
    if exist(IS_file,'file')
        mri = lab_read_data(IS_file);
    else
        mri = [];
    end
    if isfield(mri,'anatomy')
        tmp=strfind(settings.LOC_file,filesep);
        if isempty(tmp)
            if strcmp(settings.LOC_file,'data') & isfield(header,'locs')
                LOCS = header.locs;
            else
                findlocs = dir(fullfile(header.EEG_filepath,['*.' settings.LOC_file]));
                if size(findlocs,1) > 0
                    LOCS = lab_read_locs(fullfile(header.EEG_filepath,findlocs(1,1).name));
                else
                    LOCS = [];
                end
                clearvars findlocs
            end
        else
            LOCS = lab_read_locs(settings.LOC_file);
        end
        clearvars tmp
    else
        LOCS = [];
    end
    if ~isempty(LOCS) & isfield(settings,'LF')
        if settings.LF.coregmode == 3
            if ~isfield(settings.LF,'landmarks') | isempty(settings.LF.landmarks)
                settings2.LOCS = LOCS;
                settings2.indexed = mrilandmarks;
                settings2.Color = [1 1 1];
                settings2.ColorIdx = [1 0 0];
                settings2.Title = 'Select Channels for Landmarks';
                settings.LF.landmarks.landmarks = lab_plot_locs(settings2,1);
                clearvars settings2
            end
            settings2 = settings.LF.landmarks;
            settings2.forceselection = true;
            lab_mri_landmarks(mri,LOCS,settings2);
        elseif settings.coregmode == 2
            lab_coreg(IS_file,LOCS,settings.LF);
        end
    end
end

end
