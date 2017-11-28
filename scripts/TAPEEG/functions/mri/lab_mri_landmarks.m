% function to landmarks in a mri using lab_plot_orthoslides
%   Input is mri, locs (to select landmarks according to electrodes)
%
% [mri,settings] = lab_mri_landmarks(mri,LOCS,settings)
%
% written by F. Hatz 2014

function [mri,settings] = lab_mri_landmarks(mri,LOCS,settings)

disp('    Set MRI landmarks')
flagchar = 0;

if ~exist('settings','var') | ~isfield(settings,'forceselection')
    settings.forceselection = false;
end
if ~exist('LOCS','var')
    LOCS = [];
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

if isfield(mri,'landmarks')
    mrilandmarks = [];
    for i = 1:length(mri.landmarks)
        tmp = find(strcmp(LOCS.labels,mri.landmarks(i).name));
        if ~isempty(tmp)
            mrilandmarks = [mrilandmarks tmp];
        end
    end
else
    mrilandmarks = [];
end

if ~isempty(LOCS)
    if isfield(settings,'landmarks') & ~isempty(settings.landmarks)
        indexed = settings.landmarks;
    else
        settings2.LOCS = LOCS;
        settings2.indexed = mrilandmarks;
        settings2.Color = [1 1 1];
        settings2.ColorIdx = [1 0 0];
        settings2.Title = 'Select Channels for Landmarks';
        indexed = lab_plot_locs(settings2,1);
    end
end

list = {};
if isfield(mri,'landmarks') & ~isempty(mri.landmarks)
    for i = 1:length(mri.landmarks)
        list{i,1} = mri.landmarks(i).name;
    end
end

doselection = false;
landmarks = [];
for i = 1:length(indexed)
    landmarks(i).name = LOCS.labels{indexed(i)};
    tmp = find(strcmp(list,landmarks(i).name));
    if ~isempty(tmp)
        landmarks(i).pnt = mri.landmarks(tmp(1)).pnt;
        landmarks(i).mode = mri.landmarks(tmp(1)).mode;
        landmarks(i).calcfactor = mri.landmarks(tmp(1)).calcfactor;
        landmarks(i).calc1 = mri.landmarks(tmp(1)).calc1;
        landmarks(i).calc2 = mri.landmarks(tmp(1)).calc2;
    else
        landmarks(i).pnt = [];
        landmarks(i).mode = 'fixed';
        landmarks(i).calcfactor = [];
        landmarks(i).calc1 = [];
        landmarks(i).calc2 = [];
        doselection = true;
    end
end
mri.landmarks = landmarks;

if doselection == true | settings.forceselection == true
    settings2.fhandle = figure('Color',[1 1 1],'InvertHardCopy','off','NumberTitle','off', ...
        'Menubar','none','Name','Correct Orientation to RAS');
    pos = get(0,'ScreenSize');
    pos = [100 (pos(4)-800) 700 700];
    if pos(2) < 0
        pos(2) = 100;
    end
    set(settings2.fhandle,'Position',pos);
    settings2.setlandmarks = true;
    mri = lab_plot_orthoslides(mri,settings2);
    close(settings2.fhandle);
    if isfield(mri,'landmarks') & isfield(mri,'storelandmarks') & mri.storelandmarks == true
        landmarks = mri.landmarks;
        if exist('MRI_filepath','var')
            disp(['    save MRI landmarks to file ' MRI_fileS '.lmrk'])
            save(fullfile(MRI_filepath,[MRI_fileS '.lmrk']),'landmarks');
        end
    else
        mri = [];
    end
end

if flagchar == 1
    if ~isempty(mri)
        mri = MRI_file;
    else
        mri = '';
    end
end

end