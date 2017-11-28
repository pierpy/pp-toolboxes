% Set mri-orientation to RAS using lab_plot_orthoslides
%
% [mri,settings] = lab_correct_mri2RAS(mri,settings)
%
% written by F. Hatz 2014

function [mri,settings] = lab_correct_mri2RAS(mri,settings)

disp('Correct MRI-orientation to RAS')

if ~exist('settings','var')
    settings = [];
end
if ~exist('mri','var')
    [mri,MRI_file] = lab_read_mri;
    [~,MRI_filepath,~,MRI_fileS] = lab_filename(MRI_file);
elseif ischar(mri)
    flagchar = 1;
    [mri,MRI_file] = lab_read_mri(mri);
    [~,MRI_filepath,~,MRI_fileS] = lab_filename(MRI_file);
elseif ~isfield(mri,'anatomy')
    mri = [];
end
if isempty(mri)
    return
end


if exist('MRI_filepath','var') & exist(fullfile(MRI_filepath,[MRI_fileS '_RAS.hdr']),'file')
    disp('  read RAS oriented MRI from previous run')
    mri = lab_read_mri(fullfile(MRI_filepath,[MRI_fileS '_RAS.hdr']));
else
    settings.showorient = true;
    if exist('MRI_fileS','var')
        settings.fhandle = figure('Color',[1 1 1],'InvertHardCopy','off','NumberTitle','off', ...
        'Menubar','none','Name',['Correct Orientation to RAS: ' MRI_fileS]);
    else
        settings.fhandle = figure('Color',[1 1 1],'InvertHardCopy','off','NumberTitle','off', ...
            'Menubar','none','Name','Correct Orientation to RAS');
    end
    pos = get(0,'ScreenSize');
    pos = [100 (pos(4)-800) 700 700];
    if pos(2) < 0
        pos(2) = 100;
    end
    set(settings.fhandle,'Position',pos);
    flag = 0;
    while flag == 0
        [mri,settings] = lab_plot_orthoslides(mri,settings);
        answer = questdlg('Is the orientation RAS?','Orientation','Cancel','Yes','No','No');
        if strcmp(answer,'Cancel')
            mri = [];
            close(settings.fhandle);
            return
        elseif strcmp(answer,'Yes')
            if exist('MRI_filepath','var')
                lab_write_hdr(fullfile(MRI_filepath,[MRI_fileS '_RAS.hdr']),mri);
                if exist('flagchar','var')
                    mri = fullfile(MRI_filepath,[MRI_fileS '_RAS.hdr']);
                end
            end
            close(settings.fhandle);
            flag = 1;
        end
    end
end

end