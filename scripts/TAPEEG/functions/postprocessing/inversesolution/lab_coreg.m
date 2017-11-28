% Coregistration of surface points to surface
%
% [Transform,LOCS,digits] = lab_coreg(MRI_file,LOC_file,duplicates)
%
% MRI_file      =  Path to mri-file, structure with mri-data or
%                  vertices (nx3) of scalp-mesh
% LOC_file      =  xyz- or els-file with locations or fif-file
%                  (if fif-file sensor-locs and digitizer-points are
%                   processed)
% settings
%     coregmode     =  1 = automatic
%                      2 = interactive (fieldtrip code)
%                      3 = Landmarks
%     maxdist       =  maximal distance for point on nose (only for digitizer points)
%                      default = 1
%     correctdigits =  show window to remove manually outliers (only for digitizer points)
%                      default = 1
%
% written by H. Bousleiman / F. Hatz 2013

function [LOCS,Transform,settings] = lab_coreg(MRI_file,LOC_file,settings)

if ~exist('MRI_file','var')
    [settings,skipprocessing] = lab_set_coreg;
    if skipprocessing == 1
        LOCS = [];
        Transform = [];
        return
    else
        pause(0.2);
    end
    MRI_file = settings.MRI_file;
    LOC_file = settings.LOC_file;
end

if ~exist('MRI_file','var')
    [MRI_file,MRI_filepath]=uigetfile('*.hdr;*.nii;*.fif;*.fiff','Select MRI file');
    cd(MRI_filepath);
    [~,~,~,MRI_fileS] = lab_filename(MRI_file);
    warning off %#ok<WNOFF>
    mkdir(fullfile(MRI_filepath,'MriSeg'));
    warning off %#ok<WNOFF>
    MRI_file_out = fullfile(fullfile(MRI_filepath,'MriSeg'),MRI_fileS);
elseif ischar(MRI_file)
    [MRI_file,MRI_filepath,~,MRI_fileS] = lab_filename(MRI_file);
    cd(MRI_filepath);
    warning off %#ok<WNOFF>
    mkdir(fullfile(MRI_filepath,'MriSeg'));
    warning off %#ok<WNOFF>
    MRI_file_out = fullfile(fullfile(MRI_filepath,'MriSeg'),MRI_fileS);
end
if ~exist('LOC_file','var')
    [LOC_file,LOC_filepath]=uigetfile('*.xyz;*.els;*.fif;*.fiff','Select location file');
    [~,~,LOC_format,LOC_fileS] = lab_filename(LOC_file);
elseif ischar(LOC_file);
    [LOC_file,LOC_filepath,LOC_format,LOC_fileS] = lab_filename(LOC_file);
end

if ~exist('LOC_format','var') | ~strcmp(LOC_format','xyz')
    LOC_format = 'els';
end

if ischar(LOC_file)
    if exist(['lab_read_' LOC_format]) == 2
        eval(['tmp = nargout(@lab_read_' LOC_format ');']);
        if tmp > 1
            disp(['    Read locs from ' fullfile(LOC_filepath,LOC_file)])
            eval(['[~,tmp] = lab_read_' LOC_format '(fullfile(LOC_filepath,LOC_file));']);
            if isfield(tmp,'locs')
                LOCS = tmp.locs;
            else
                disp('Abort: no valid LOC-file')
                LOCS = [];
                Transform = [];
                return
            end
        else
            LOCS = lab_read_locs(fullfile(LOC_filepath,LOC_file));
        end
        clearvars tmp
    else
        disp('Abort: no valid location-info')
        LOCS = [];
        Transform = [];
        return
    end
    clearvars LOCformat
elseif isnumeric(LOC_file) & size(LOC_file,2) == 3
    LOCS.chanpos = LOC_file;
elseif isstruct(LOC_file) & (isfield(LOC_file,'x') | isfield(LOC_file,'chanpos'))
    LOCS = LOC_file;
else
    disp('Abort: no valid location-info')
    LOCS = [];
    Transform = [];
    return
end

if isfield(LOCS,'x')
    if isfield(LOCS,'chanpos')
        LOCS = rmfield(LOCS,'chanpos');
    end
    LOCS.chanpos(:,1) = LOCS.x';
    LOCS.chanpos(:,2) = LOCS.y';
    LOCS.chanpos(:,3) = LOCS.z';
end
if ~isfield(LOCS,'label')
    LOCS.label = cellstr(num2str((1:size(LOCS.chanpos,1))'));
    LOCS.unit = 'mm';
    LOCS.elecpos = LOCS.chanpos;
end

if isfield(LOCS,'x')
    numlocs = length(LOCS.x);
elseif isfield(LOCS,'chanpos')
    numlocs = size(LOCS.chanpos,1);
end
if isfield(settings,'deletelocs') & ~isempty(settings.deletelocs)
    numlocs = length(setdiff(1:numlocs,settings.deletelocs));
end

if exist('MRI_file_out','var')
    cd(MRI_filepath);
    warning off %#ok<WNOFF>
    mkdir(fullfile(MRI_filepath,'Coreg'));
    warning off %#ok<WNOFF>
    LOC_file_out = fullfile(fullfile(MRI_filepath,'Coreg'),[MRI_fileS '_E' num2str(numlocs)]);
elseif ischar(LOC_file)
    cd(LOC_filepath);
    warning off %#ok<WNOFF>
    mkdir(fullfile(LOC_filepath,'Coreg'));
    warning off %#ok<WNOFF>
    LOC_file_out = fullfile(fullfile(LOC_filepath,'Coreg'),[LOC_fileS '_Coreg']);
end

if exist('LOC_file_out','var') & exist([LOC_file_out '.xyz'],'file')
    disp('    Read coregistered locs from previous run')
    LOCS = lab_read_locs([LOC_file_out '.xyz']);
    LOCS.chanpos(:,1) = LOCS.x';
    LOCS.chanpos(:,2) = LOCS.y';
    LOCS.chanpos(:,3) = LOCS.z';
    LOCS.elecpos = LOCS.chanpos;
    LOCS.label = LOCS.labels';
    LOCS.unit = 'mm';
    Transform = [];
    return
end

if ~exist('settings','var') | ~isfield(settings,'maxdist')
    maxdist = 1;
else
    maxdist = settings.maxdist;
end
if ~exist('settings','var') | ~isfield(settings,'correctdigits')
    correctdigits = 1;
else
    correctdigits = settings.correctdigits;
end
if ~exist('settings','var') | ~isfield(settings,'coregmode')
    coregmode = 1;
else
    coregmode = settings.coregmode;
end

if ischar(MRI_file)
    if ~exist([MRI_file_out '_scalp.hdr'],'file')
        cfgtmp = [];
        cfgtmp.mrifile = fullfile(MRI_filepath,MRI_file);
        cfgtmp.SEGcorrect = 1;
        cfgtmp.SEGprobmaps = 0;
        mri = lab_segment_mri(cfgtmp);
        clearvars cfgtmp
    else
        disp(['    Read Scalp-Info from ' MRI_file_out '_scalp.hdr'])
        mri = lab_read_mri(fullfile(MRI_filepath,MRI_file));
        mritmp = lab_read_mri([MRI_file_out '_scalp.hdr']);
        mri.scalp = logical(mritmp.anatomy);
        clearvars mritmp
    end
elseif isstruct(MRI_file) & isfield(MRI_file,'scalp')
    mri = MRI_file;
    mri = lab_set_mricoord(mri);
    MRI_file = [];
end
if exist('mri','var') & isstruct(mri) & isfield(mri,'scalp')
    cfgtmp = [];
    cfgtmp.nvert = 20000;
    [bnd,cfgtmp] = lab_mesh_segm(mri,cfgtmp);
    if isfield(cfgtmp,'tissue')
        tmp = find(strcmp(cfgtmp.tissue,'scalp'));
        if isempty(tmp)
            tmp = 1;
        end
    else
        tmp = 1;
    end
    [scalp,scalptri] = tess_refine(bnd(tmp).pnt,bnd(tmp).tri,5);
    [scalp,~] = tess_refine(scalp,scalptri,5);
    clearvars cfgtmp
elseif isnumeric(MRI_file) & (size(MRI_file,2) == 3 | size(MRI_file,1) == 3)
    scalp = MRI_file;
    if size(MRI_file,1) == 3
        scalp = scalp';
    end
    MRI_file = [];
else
    disp('Abort: no valid mri-info')
    return
end

if isfield(LOCS,'digits')
    digits = LOCS.digits;
else
    digits = [];
end

if coregmode == 3
    if ~isempty(digits)
        disp('    Landmarks coregistration for MEG data with digits not possible, change to automatic mode')
        coregmode = 1;
    elseif ~exist('mri','var') | ~isstruct(mri) | ~isfield(mri,'anatomy')
        disp('    Landmarks coregistration without MRI not possible, change to automatic mode')
        coregmode = 1;
    end
end

if coregmode == 1
    if ~isempty(digits)
        % show window to delete outliers
        if correctdigits == 1
            digits = lab_correct_digits(digits);
            pause(0.2);
            if exist('LOC_filepath','var')
                LOCS.digits = digits;
                lab_write_locs(fullfile(LOC_filepath,[LOC_fileS '.' LOC_format]),LOCS);
            end
        end
        
        % find digits of nose (lowest 16% of digits)
        tmp = max(digits(:,3)) - min(digits(:,3));
        tmp = min(digits(:,3)) + tmp/6;
        digitsNose = digits(digits(:,3)<=tmp,:);
        clearvars tmp
        
        % delete scalp 1cm below nose
        tmp = scalp(find(scalp(:,2)==max(scalp(:,2)),1,'first'),:);
        scalp = scalp(scalp(:,3)>=tmp(1,3)-10,:);
        clearvars tmp
        
        disp(['    Coregister using icp (with maximal distance ' num2str(maxdist) ' for nose)'])
        dist = maxdist + 1;
        counter = 1;
        digitstmp = digits;
        while dist > maxdist
            if counter > 40
                disp('Abort, coregistration not possible')
                LOCS= [];
                Transform = [];
                return
            end
            % do icp
            disp(['     do icp with ' num2str(counter) ' nose replicates'])
            [R, T] = icp(scalp', digitstmp',1000,'Matching','kDtree','WorstRejection',0.01);
            Transform = eye(4);
            Transform(1:3,1:3) = R;
            Transform(1:3,4) = T;
            
            % check noise
            tmp = (Transform * [digitsNose ones(size(digitsNose,1),1)]')';
            [~,dist] = knnsearch(scalp,tmp(:,1:3));
            dist = sort(dist);
            dist = min(dist(1:round(length(dist)/2)));
            clearvars tmp
            disp(['     distance to nose: ' num2str(dist)])
            
            % add some more noise points
            digitstmp = [digitstmp;digitsNose];
            counter = counter + 1;
        end
        clearvars dist digitstmp
    else
        disp('    Coregister using icp')
        % do icp
        [R, T] = icp(scalp',LOCS.chanpos',1000,'Matching','kDtree');
        Transform = eye(4);
        Transform(1:3,1:3) = R;
        Transform(1:3,4) = T;
    end
elseif coregmode == 2
    disp('    Coregister using interactive mode (fieldtrip)')
    cfg.headshape = scalp;
    cfg.method = 'interactive';
    if ~isempty(digits)
        cfg.elec.chanpos = digits;
        cfg.elec.label = cellstr(num2str((1:size(digits,1))'));
        cfg.elec.unit = 'mm';
        cfg.elec.elecpos = cfg.elec.chanpos;
    else
        cfg.elec = LOCS;
    end
    warning off %#ok<WNOFF>
    % spmpath = spm('dir');
    % cfg.template = [spmpath,filesep,'templates',filesep,'T1.nii'];
    % cfg.tpm   = {fullfile(spmpath,'tpm','grey.nii'),...
    %     fullfile(spmpath,'tpm','white.nii'),fullfile(spmpath,'tpm','csf.nii')};
    norm = ft_electroderealign(cfg);
    pause(0.2);
    warning on %#ok<WNON>
    Transform = norm.m;
elseif coregmode == 3
    if ischar(MRI_file)
        scalptmp = mri.scalp;
        if isfield(settings,'landmarks') & isfield(settings.landmarks,'mrilandmarks') & ...
                exist(settings.landmarks.mrilandmarks,'file')
            mrifiletmp = lab_match_mrilandmarks2mri(fullfile(MRI_filepath,MRI_file),settings.landmarks.mrilandmarks);
            mri = lab_read_mri(mrifiletmp);
        end
        mri.scalp = scalptmp;
        clearvars mrifiletmp scalptmp
    else
        if isfield(settings,'landmarks') & isfield(settings.landmarks,'mrilandmarks') && ...
                exist(settings.landmarks.mrilandmarks,'file')
            mri = lab_match_mrilandmarks2mri(mri,settings.landmarks.mrilandmarks);
        end
    end
    if isfield(LOCS,'chanpos') & ~isfield(LOCS,'x')
        LOCS.x = LOCS.chanpos(:,1)';
        LOCS.y = LOCS.chanpos(:,2)';
        LOCS.z = LOCS.chanpos(:,3)';
    end
    if isfield(LOCS,'label') & ~isfield(LOCS,'labels')
        LOCS.labels = LOCS.label(:)';
    end
    [mri,LOCS,settings.landmarks] = lab_warp_locs2mri(mri,LOCS,scalp,settings);
    if isfield(mri,'landmarks') & exist('MRI_filepath','var')
        landmarks = mri.landmarks;
        save(fullfile(MRI_filepath,[MRI_fileS '.lmrk']),'landmarks');
    end
    Transform = [];
    LOCS  = lab_locs2sph(LOCS);
    if isfield(LOCS,'chanpos')
        LOCS.chanpos = [LOCS.x' LOCS.y' LOCS.z'];
    end
end

if coregmode < 3
    % apply transfrom to LOCS and digits
    if ~isempty(digits)
        digits = (Transform * [digits ones(size(digits,1),1)]')';
        digits = digits(:,1:3);
        LOCS.digits = digits;
    end
    if isfield(LOCS,'chanpos')
        chans =  (Transform * [LOCS.chanpos ones(size(LOCS.chanpos,1),1)]')';
        LOCS.chanpos = chans(:,1:3);
    end
    if isfield(LOCS,'grad')
        if isfield(LOCS.grad,'unit') & ~strcmp(LOCS.grad,'mm')
            LOCS.grad = ft_convert_units(LOCS.grad,'mm');
        end
        LOCS.grad.coilori = (Transform(1:3,1:3) * LOCS.grad.coilori')';
        coilpos =  (Transform * [LOCS.grad.coilpos ones(size(LOCS.grad.coilpos,1),1)]')';
        LOCS.grad.coilpos = coilpos(:,1:3);
    end
    if isfield(LOCS,'x')
        LOCS.x = LOCS.chanpos(:,1)';
        LOCS.y = LOCS.chanpos(:,2)';
        LOCS.z = LOCS.chanpos(:,3)';
        LOCS  = lab_locs2sph(LOCS);
    end
end

% Delete selected electrodes, if necessary
if isfield(settings,'deletelocs') & ~isempty(settings.deletelocs)
    disp('    Delete electrode positions after coregistration')
    if isfield(LOCS,'x')
        numlocs = length(LOCS.x);
    elseif isfield(LOCS,'chanpos')
        numlocs = size(LOCS.chanpos,1);
    end
    selection = setdiff(1:numlocs,settings.deletelocs);
    if isfield(LOCS,'x')
        LOCS.x = LOCS.x(1,selection);
        LOCS.y = LOCS.y(1,selection);
        LOCS.z = LOCS.z(1,selection);
        LOCS.labels = LOCS.labels(1,selection);
        LOCS  = lab_locs2sph(LOCS);
    elseif isfield(LOCS,'chanpos')
        LOCS.chanpos = LOCS.chanpos(selection,:);
    end
end

% write results
if exist('LOC_file_out','var')
    Invisible = 1;
    if exist('Transform','var') & ~isempty(Transform)
        dlmwrite([LOC_file_out '_Tlocs2scalp.txt'],Transform,'delimiter','\t','precision',6,'newline','pc');
        dlmwrite([LOC_file_out '_Tscalp2locs.txt'],inv(Transform),'delimiter','\t','precision',6,'newline','pc');
    end
    lab_write_locs([LOC_file_out '.xyz'],LOCS);
    lab_write_spi([LOC_file_out '_scalp.pts'],scalp,1);
else
    Invisible = 0;
end
% if exist('LOC_filepath','var')
%     lab_write_locs(fullfile(LOC_filepath,[LOC_fileS '.' LOC_format]),LOCS);
% end
if exist('mri','var') & isstruct(mri) & isfield(mri,'anatomy') & exist('LOC_file_out','var') & ...
        exist('Transform','var') & ~isempty(Transform)
    mri.transform = inv(Transform) * mri.transform;
    ft_write_mri([LOC_file_out '_coreg.nii'],mri.anatomy,'dataformat','nifti','transform',mri.transform);
    ElektaLandmarks = (Transform * [-70 0 0 1;0 100 0 1;70 0 0 1]')';
    ElektaLandmarks = ElektaLandmarks(:,1:3);
    dlmwrite([LOC_file_out '_Landmarks.txt'],ElektaLandmarks,'delimiter','\t','precision',4,'newline','pc');
end

% plot results
plot.facecolor = [0.75 0.75 0.75];
plot.alpha = 0.7;
plot.invisible = 1;
plot.plotfaces = true;
plot.plotedges = false;
f = lab_plot_mesh(scalp,plot);
if ~isempty(digits)
    plot.dotscolor = [1 0 0];
    plot.plotdots = true;
    lab_plot_mesh(digits,plot);
end
if ~isempty(LOCS)
    if ~isempty(digits) & isfield(LOCS,'grad')
        plot.plotdots = false;
        plot.facecolor = [0 0 1];
        lab_plot_mesh(LOCS.grad,plot);
    else
        plot.dotscolor = [0 0 1];
        plot.plotdots = true;
        lab_plot_mesh(LOCS,plot);
    end
end
if exist('LOC_file_out','var')
    warning off %#ok<WNOFF>
    lab_print_figure([LOC_file_out '_coregS.jpg'],f);
    view(gca,[0 90]);
    lab_print_figure([LOC_file_out '_coregT.jpg'],f);
    view(gca,[180 0]);
    lab_print_figure([LOC_file_out '_coregF.jpg'],f);
    close(f);
    warning on %#ok<WNON>
end

% Calculate sphere origin for Elekta beamformer
if exist('LOC_file_out','var') & isfield(settings,'findsphere') & ...
        settings.findsphere == 1 & settings.coregmode ~= 3
    disp('    Calculate sphere origin for elekta beamformer')
    cfgtmp.method = 'singlesphere';
    bnd.pnt = scalp;
    bnd.tri = lab_pnt2faces(bnd.pnt);
    voltmp = ft_prepare_headmodel(cfgtmp,bnd);
    sphere.origin = (inv(Transform) * [voltmp.o 1]')';
    sphere.origin = sphere.origin(1,1:3);
    sphere.origin(1) = -sphere.origin(1);
    sphere.radius= voltmp.r;
    clearvars voltmp cfgtmp
    dlmwrite([LOC_file_out '_Sphere.txt'],sphere.origin,'delimiter','\t','precision',6,'newline','pc');
end
