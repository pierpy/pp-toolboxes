% Segment mri to scalp, skull and brain (optionally csf, gray and white)
%
% mri = lab_segment_mri(settings)
%   settings.mrifile = MRI_file
%
% written by F. Hatz 2012

function mri = lab_segment_mri(settings)

disp('    Segment mri using fieldtrip/spm code')

if exist('settings','var')
    if ischar(settings)
        tmp = settings;
        clearvars settings
        if exist(tmp,'file')
            settings.mrifile = tmp;
        end
        clearvars tmp
    elseif isstruct(settings) & isfield(settings,'anatomy')
        mri = settings;
        if ~isfield(mri,'coordsys')
            mri.coordsys = 'ras';
        end
        clearvars settings
        settings.mrifile = [];
    end
end
if ~exist('settings','var') | ~isfield(settings,'mrifile')
    [filename,filepath] = uigetfile('*.hdr;*.fif;*.nii','select MRI-file');
    settings.mrifile = fullfile(filepath,filename);
    clearvars filename filepath
end
if ~isempty(settings.mrifile) & ischar(settings.mrifile)
    [~,mripath,~,mrifileS] = lab_filename(settings.mrifile);
    warning off
    mkdir(fullfile(mripath,'MriSeg'));
    warning on
    mrifileS = fullfile(fullfile(mripath,'MriSeg'),mrifileS);
    clearvars mripath
else
    mrifileS = [];
end

if ~isfield(settings,'SEGcorrect')
    settings.SEGcorrect = 1;
end
if ~isfield(settings,'SEGprobmaps')
    settings.SEGprobmaps = 1;
end
if ~isfield(settings,'SEGbrain')
    settings.SEGbrain = 1;
end

if ~exist('mri','var')
    if ischar(settings.mrifile);
        [mri,settings.mrifile] = lab_read_mri(settings.mrifile);
    elseif isstruct(settings.mrifile)
        mri = settings.mrifile;
    end
end

if settings.SEGbrain == 1
    if isempty(mrifileS) | ~exist([mrifileS '_scalp.hdr'],'file') |  ~exist([mrifileS '_brain.hdr'],'file') |  ...
            ~exist([mrifileS '_skull.hdr'],'file') | ~exist([mrifileS '_csf.hdr'],'file') | ...
            ~exist([mrifileS '_gray.hdr'],'file') |  ~exist([mrifileS '_white.hdr'],'file')
        cfg = [];
        cfg.output = {'gray','white','csf','skull','scalp'};
        if isfield(settings,'scalpthreshold')
            cfg.scalpthreshold = settings.scalpthreshold;
        end
        if isfield(settings,'brainthreshold')
            cfg.brainthreshold = settings.brainthreshold;
        end
        spmpath = spm('dir');
        cfg.template = [spmpath,filesep,'templates',filesep,'T1.nii'];
        cfg.tpm   = {fullfile(spmpath,'tpm','grey.nii'),...
            fullfile(spmpath,'tpm','white.nii'),fullfile(spmpath,'tpm','csf.nii')};
        warning off %#ok<WNOFF>
        segmentedmri = ft_volumesegment(cfg,mri);
        warning on %#ok<WNON>
        mri.scalp = segmentedmri.scalp;
        mri.skull = segmentedmri.skull;
        mri.csf = segmentedmri.csf;
        mri.gray = segmentedmri.gray;
        mri.white = segmentedmri.white;
        mri.brain = mri.csf;
        mri.brain(mri.gray == 1) = 1;
        mri.brain(mri.white == 1) = 1;
        clearvars segmentedmri
        % correct for errors
        if settings.SEGcorrect == 1
            mri = lab_correct_segm(mri);
        end
        if ~isempty(mrifileS)
            if ~isfield(mri,'dimensions')
                mri = lab_set_mricoord(mri);
            end
            mriout = mri;
            mriout.anatomy = uint8(mri.scalp);
            lab_write_hdr([mrifileS '_scalp.hdr'],mriout);
            mriout = mri;
            mriout.anatomy = uint8(mri.skull);
            lab_write_hdr([mrifileS '_skull.hdr'],mriout);
            mriout = mri;
            mriout.anatomy = uint8(mri.brain);
            lab_write_hdr([mrifileS '_brain.hdr'],mriout);
            mriout = mri;
            mriout.anatomy = uint8(mri.csf);
            lab_write_hdr([mrifileS '_csf.hdr'],mriout);
            mriout = mri;
            mriout.anatomy = uint8(mri.gray);
            lab_write_hdr([mrifileS '_gray.hdr'],mriout);
            mriout = mri;
            mriout.anatomy = uint8(mri.white);
            lab_write_hdr([mrifileS '_white.hdr'],mriout);
            combine = zeros(size(mri.scalp));
            combine(mri.scalp==1) = 250;
            combine(mri.skull==1) = 200;
            combine(mri.csf==1) = 150;
            combine(mri.gray==1) = 100;
            combine(mri.white==1) = 50;
            mriout.anatomy = uint8(combine);
            lab_write_hdr([mrifileS '_segments.hdr'],mriout);
        end
    else
        disp('     load previously calculated files (scalp/skull/brain)')
        tmp = lab_read_mri([mrifileS '_scalp.hdr'],1);
        mri.scalp = logical(tmp.anatomy);
        tmp = lab_read_mri([mrifileS '_skull.hdr'],1);
        mri.skull = logical(tmp.anatomy);
        tmp = lab_read_mri([mrifileS '_brain.hdr'],1);
        mri.brain = logical(tmp.anatomy);
        tmp = lab_read_mri([mrifileS '_csf.hdr'],1);
        mri.csf = logical(tmp.anatomy);
        tmp = lab_read_mri([mrifileS '_gray.hdr'],1);
        mri.gray = logical(tmp.anatomy);
        tmp = lab_read_mri([mrifileS '_white.hdr'],1);
        mri.white = logical(tmp.anatomy);
        clearvars tmp
        if ~exist([mrifileS '_segments.hdr'],'file')
            combine = zeros(size(mri.scalp));
            combine(mri.scalp==1) = 250;
            combine(mri.skull==1) = 200;
            combine(mri.csf==1) = 150;
            combine(mri.gray==1) = 100;
            combine(mri.white==1) = 50;
            mriout = mri;
            mriout.anatomy = uint8(combine);
            lab_write_hdr([mrifileS '_segments.hdr'],mriout);
        end
    end
else
    if isempty(mrifileS) | ~exist([mrifileS '_scalp.hdr'],'file') |  ...
            ~exist([mrifileS '_skull.hdr'],'file') | ~exist([mrifileS '_brain.hdr'],'file')
        cfg = [];
        cfg.output = {'brain','skull','scalp'};
        if isfield(settings,'scalpthreshold')
            cfg.scalpthreshold = settings.scalpthreshold;
        end
        if isfield(settings,'brainthreshold')
            cfg.brainthreshold = settings.brainthreshold;
        end
        spmpath = spm('dir');
        cfg.template = [spmpath,filesep,'templates',filesep,'T1.nii'];
        cfg.tpm   = {fullfile(spmpath,'tpm','grey.nii'),...
            fullfile(spmpath,'tpm','white.nii'),fullfile(spmpath,'tpm','csf.nii')};
        warning off %#ok<WNOFF>
        segmentedmri = ft_volumesegment(cfg, mri);
        warning on %#ok<WNON>
        mri.brain = segmentedmri.brain;
        mri.scalp = segmentedmri.scalp;
        mri.skull = segmentedmri.skull;
        clearvars segmentedmri
        % correct for errors
        if settings.SEGcorrect == 1
            mri = lab_correct_segm(mri);
        end
        if ~isempty(mrifileS)
            if ~isfield(mri,'dimensions')
                mri = lab_set_mricoord(mri);
            end
            mriout = mri;
            mriout.anatomy = uint8(mri.scalp);
            lab_write_hdr([mrifileS '_scalp.hdr'],mriout);
            mriout = mri;
            mriout.anatomy = uint8(mri.skull);
            lab_write_hdr([mrifileS '_skull.hdr'],mriout);
            mriout = mri;
            mriout.anatomy = uint8(mri.brain);
            lab_write_hdr([mrifileS '_brain.hdr'],mriout);
            combine = zeros(size(mri.scalp));
            combine(mri.scalp==1) = 250;
            combine(mri.skull==1) = 150;
            combine(mri.brain==1) = 50;
            mriout.anatomy = uint8(combine);
            lab_write_hdr([mrifileS '_seg.hdr'],mriout);
        end
    else
        disp('     load previously calculated files (scalp/skull/brain)')
        tmp = lab_read_mri([mrifileS '_scalp.hdr'],1);
        mri.scalp = logical(tmp.anatomy);
        tmp = lab_read_mri([mrifileS '_skull.hdr'],1);
        mri.skull = logical(tmp.anatomy);
        tmp = lab_read_mri([mrifileS '_brain.hdr'],1);
        mri.brain = logical(tmp.anatomy);
        clearvars tmp
        if ~exist([mrifileS '_seg.hdr'],'file')
            combine = zeros(size(mri.scalp));
            combine(mri.scalp==1) = 250;
            combine(mri.skull==1) = 150;
            combine(mri.brain==1) = 50;
            mriout = mri;
            mriout.anatomy = uint8(combine);
            lab_write_hdr([mrifileS '_segments.hdr'],mriout);
        end
    end
end

if settings.SEGprobmaps == 1
    if (isempty(mrifileS) | ~exist([mrifileS '_csf_tpm.hdr'],'file') |  ...
            ~exist([mrifileS '_gray_tpm.hdr'],'file') | ~exist([mrifileS '_white_tpm.hdr'],'file'))
        cfg = [];
        spmpath = spm('dir');
        cfg.template = [spmpath,filesep,'templates',filesep,'T1.nii'];
        cfg.tpm   = {fullfile(spmpath,'tpm','grey.nii'),...
            fullfile(spmpath,'tpm','white.nii'),fullfile(spmpath,'tpm','csf.nii')};
        warning off %#ok<WNOFF>
        segmentedmri = ft_volumesegment(cfg,mri);
        warning on %#ok<WNON>
        mri.csf = segmentedmri.csf;
        mri.gray = segmentedmri.gray;
        mri.white = segmentedmri.white;
        clearvars segmentedmri
        
        % correct volumes (match with brain-volume
        if settings.SEGcorrect == 1 & isfield(mri,'brain')
            mri = lab_correct_segm(mri);
        end
        % write results
        if ~isempty(mrifileS)
            mriout = mri;
            mriout.anatomy = mri.csf;
            lab_write_hdr([mrifileS '_csf_tpm.hdr'],mriout);
            mriout = mri;
            mriout.anatomy = mri.gray;
            lab_write_hdr([mrifileS '_gray_tpm.hdr'],mriout);
            mriout = mri;
            mriout.anatomy = mri.white;
            lab_write_hdr([mrifileS '_white_tpm.hdr'],mriout);
        end
    else
        disp('     load previously calculated files (csf(gray/white)')
        tmp = lab_read_mri([mrifileS '_csf_tpm.hdr'],1);
        mri.csf = tmp.anatomy;
        tmp = lab_read_mri([mrifileS '_gray_tpm.hdr'],1);
        mri.gray = tmp.anatomy;
        tmp = lab_read_mri([mrifileS '_white_tpm.hdr'],1);
        mri.white = tmp.anatomy;
        clearvars tmp
    end
end

return

