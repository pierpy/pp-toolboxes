% Processing of MRI's, script called by lab_tapeeg
%
% [cfg,Result] = lab_process_mri(cfg)
%
% Written by F. Hatz 2013 Neurology Basel

function [cfg,Result] = lab_process_mri(cfg)

Result = [];
if ~exist('cfg','var') | ~isfield(cfg,'EEG_file')
    [cfg.EEG_file,cfg.EEG_filepath] = uigetfile('*.hdr;*.nii;*.fif','Select mri');
    if cfg.EEG_file == 0
        return
    end
end
[~,~,~,cfg.Output_fileS] = lab_filename(cfg.EEG_file);
cfg.Output_filepath = cfg.EEG_filepath;

if ~isfield(cfg,'MRI')
    [cfg,skipprocessing] = lab_set_process_mri(cfg);
    if skipprocessing == 1
        return
    end
end
cfg.MRI.mrifile = fullfile(cfg.EEG_filepath,cfg.EEG_file);

% convert MRI to '.hdr'
[~,~,Mformat] = lab_filename(cfg.MRI.mrifile);
if ~strcmp(Mformat,'hdr')
    [~,cfg.MRI.mrifile] = lab_read_mri(cfg.MRI.mrifile);
end

% correct orientation
if isfield(cfg.MRI,'orientation') & length(cfg.MRI.orientation) == 3
    [cfg.MRI.mrifile,cfg.MRI.orientation] = lab_correct_mri_orient(cfg.MRI.mrifile,cfg.MRI.orientation);
end

% correct contrast
if isfield(cfg.MRI,'CONTRAST') & ~isempty(cfg.MRI.CONTRAST)
    [cfg.MRI.mrifile,settings] = lab_adjust_mri_contrast(cfg.MRI.mrifile,cfg.MRI.CONTRAST);
    cfg.MRI.CONTRAST = settings;
    clearvars settings
end

% anisotropic filtering
if isfield(cfg.MRI,'FILT') & ~isempty(cfg.MRI.FILT)
    [cfg.MRI.mrifile,settings] = lab_filter_mri(cfg.MRI.mrifile,cfg.MRI.FILT);
    cfg.MRI.FILT = settings;
    clearvars settings
end

% match mri to template
if ~isempty(cfg.MRI.template)
    [cfg.MRI.mrifile] = lab_match_mri(cfg.MRI.mrifile,cfg.MRI.template,cfg.MRI.matchmode);
end

% segment mri
if cfg.MRI.SEGbrain == 1 | cfg.MRI.SEGprobmaps == 1
    mri = lab_segment_mri(cfg.MRI);
end

% calculate leadfield
if isfield(cfg.MRI,'LF') & ~isempty(cfg.MRI.LF) & exist(cfg.MRI.LOC_file,'file') & ~isempty(cfg.MRI.SPI)
    locs = lab_read_locs(cfg.MRI.LOC_file);
    [slocs,cfg.MRI.SPI] = lab_create_sp(cfg.MRI.mrifile,cfg.MRI.SPI);
    lab_compute_leadfield(cfg.MRI.LF,slocs,locs);
end

% collect voxels
if ~isempty(cfg.MRI.collectsegm) & isfield(mri,cfg.MRI.collectsegm) & cfg.MRI.collectthreshold > 0
    if ~isempty(cfg.MRI.atlas)
        if ~isfield(cfg.MRI,'AtlasTemplate') 
            cfg.MRI.AtlasTemplate = '';
        end
        if ~isempty(cfg.MRI.template)
            atlas = lab_match_atlas2mri(cfg.MRI.atlas,cfg.MRI.template,cfg.MRI.AtlasTemplate);
        else
            atlas = lab_match_atlas2mri(cfg.MRI.atlas,cfg.MRI.mrifile,cfg.MRI.AtlasTemplate);
        end
    else
        atlas = [];
    end
    collect = mri.(cfg.MRI.collectsegm);
    collect = collect(collect>=cfg.MRI.collectthreshold);
    Result.Segment = cfg.MRI.collectsegm;
    Result.TotalCollect = sum(collect(:));
    xlsout = {['Total ' cfg.MRI.collectsegm],Result.TotalCollect};
    if ~isempty(atlas)
        regions = unique(atlas);
        for i = regions
            tmp = collect(atlas.anatomy==i);
            Result.Regions(i) = sum(tmp(:));
            xlsout = cat(1,xlsout,{['Region' num2str(i,'%03d')],Result.Regions(i)});
        end
    end
    xlsoutfile = fullfile(cfg.Output_filepath,[cfg.Output_fileS '_' cfg.MRI.collectsegm '.xls']);
    lab_write_xls(xlsoutfile,xlsout);
end

return