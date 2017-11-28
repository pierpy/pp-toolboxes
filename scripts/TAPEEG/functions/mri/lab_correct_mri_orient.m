% Correct mri-orientation to RAS, input is mri-structure and orientation of
% input mri in format [x y z] ( RAS = [1 2 3] )
%
% [mri,orient] = lab_correct_mri_orient(mri,orient)

function [mri,orient] = lab_correct_mri_orient(mri,orient)

disp('    Correct orientation of MRI')

if ~exist('mri','var') | isempty(mri)
    orient = [];
    mri = [];
    return
end

if ~exist('orient','var') | isempty(orient)
    [orient,skipprocessing] = lab_get_orient;
    if skipprocessing == 1
        return
    end
end

if orient(1) == 1 & orient(2) == 2 & orient(3) == 3
    disp('     MRI already in RAS')
    return
end

if ischar(mri)
    [mri,mri_file] = lab_read_mri(mri);
    if isempty(mri)
        disp('     Abort: no valid MRI-file')
        return
    end
    vol = mri.anatomy;
    isnifti = false;
elseif isfield(mri,'anatomy')
    vol = mri.anatomy;
    isnifti = false;
elseif isfield(mri,'img')
    vol = mri.img;
    isnifti = true;
else
    return
end

%  calculate after flip orient
rot_orient = mod(orient + 2, 3) + 1;

%  do flip:
flip_orient = orient - rot_orient;
for i = 1:3
    if flip_orient(i)
        if exist('flipdim') == 2 %#ok<EXIST>
            vol = flipdim(vol,i);
        elseif exist('flip') == 2 %#ok<EXIST>
            vol = flip(vol,i);
        end
    end
end

%  get index of orient (do inverse)
[~,rot_orient] = sort(rot_orient);

%  do rotation:
vol = permute(vol,[rot_orient 4 5 6]);

if isnifti == true
    mri.img = vol;
    
    %  rotate resolution, or 'dim'
    new_dim = mri.hdr.dime.dim(2:4);
    new_dim = new_dim(rot_orient);
    mri.hdr.dime.dim(2:4) = new_dim;
    
    %  rotate voxel_size, or 'pixdim'
    dimension = mri.hdr.dime.pixdim(2:4);
    dimension = dimension(rot_orient);
    mri.hdr.dime.pixdim(2:4) = dimension;
    
    %  re-calculate originator
    originator = mri.hdr.hist.originator([1:3]);
    originator = originator(rot_orient);
    flip_orient = flip_orient(rot_orient);
    for i = 1:3
        if flip_orient(i) & ~isequal(double(originator(i)), 0)
            originator(i) = int16(double(new_dim(i)) - double(originator(i)) + 1);
        end
    end
    mri.hdr.hist.originator(1:3) = originator;
    
    % Set description
    if length(mri.hdr.hist.descrip) > 3
        mri.hdr.hist.descrip(end-3:end) = '_RAS';
    else
        mri.hdr.hist.descrip = [mri2.hdr.hist.descrip '_RAS'];
    end
else
    mri.anatomy = vol;
    
    mri.dim = mri.dim(1,rot_orient);
    mri.dimensions = mri.dim;
    if isfield(mri,'coordsys')
        mri.coordold = mri.coordsys;
    end
    mri.coordsys = 'ras';
    if isfield(mri,'originator')
        originator = mri.originator;
    else
        originator = mri.transform(1:3,4)';
    end
    originator = originator(rot_orient);
    flip_orient = flip_orient(rot_orient);
    for i = 1:3
        if flip_orient(i) & ~isequal(double(originator(i)), 0)
            originator(i) = int16(double(mri.dim(i)) - double(originator(i)) + 1);
        end
    end
    mri.originator = originator;
    transform = mri.transform;
    transform(1:3,1:3) = mri.transform(rot_orient,rot_orient);
    transform(1:3,4) = originator';
    mri.transform = transform;
end

if exist('mri_file','var')
    [~,mri_filepath,~,mri_fileS] = lab_filename(mri_file);
    lab_write_hdr(fullfile(mri_filepath,[mri_fileS '_RAS.hdr']),mri);
    mri = fullfile(mri_filepath,[mri_fileS '_RAS.hdr']);
end

end