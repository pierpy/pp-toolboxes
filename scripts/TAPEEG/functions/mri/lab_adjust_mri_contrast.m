function [mri,settings] = lab_adjust_mri_contrast(mri,settings,dostore)

disp('    Set contrast of MRI-file')

if ~exist('dostore','var')
    dostore = true;
end
if ~exist('settings','var')
    settings = [];
end
if ischar(mri)
    [mri,mri_file] = lab_read_mri(mri);
    vol = mri.anatomy;
elseif ~exist('mri','var')
    [mri,mri_file] = lab_read_mri;
    vol = mri.anatomy;
elseif isfield(mri,'anatomy')
    vol = mri.anatomy;
elseif isfield(mri,'img')
    vol = mri.img;
else
    mri = [];
    return
end

if isfield(settings,'minval') & ~isempty(settings.minval) & isfield(settings,'maxval') & ~isempty(settings.maxval)
    disp(['     set Min Value to ' num2str(settings.minval) ' & Max Value to ' num2str(settings.maxval)])
    vol(vol>settings.maxval) = settings.maxval;
    vol = vol - settings.minval;
    vol(vol<0) = 0;
    if isfield(mri,'img')
        mri.img = vol;
    else
        mri.anatomy = vol;
    end
elseif isfield(mri,'anatomy')
    settings2.fhandle = figure('Color',[1 1 1],'InvertHardCopy','off','NumberTitle','off', ...
        'Menubar','none','Name','Correct Orientation to RAS');
    pos = get(0,'ScreenSize');
    pos = [100 (pos(4)-800) 700 700];
    if pos(2) < 0
        pos(2) = 100;
    end
    set(settings2.fhandle,'Position',pos);
    settings2.setmaxmin = true;
    mri = lab_plot_orthoslides(mri,settings2);
    close(settings2.fhandle);
    if isfield(mri,'minval') & isfield(mri,'maxval')
        settings.minval = mri.minval;
        settings.maxval = mri.maxval;
    end
else
    disp('    Abort: selection of minimal and maximal value with input mir not possible')
    mri = [];
end

if exist('mri_file','var') & dostore == true
    [~,mri_filepath,~,mri_fileS] = lab_filename(mri_file);
    lab_write_hdr(fullfile(mri_filepath,[mri_fileS '_Correct.hdr']),mri);
    mri = fullfile(mri_filepath,[mri_fileS '_Correct.hdr']);
end

end