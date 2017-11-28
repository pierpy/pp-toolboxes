% Helper script for lab_stitching
%
% written by F. Hatz 2012

function cfg = lab_prepare_stitching_auto(cfg,header)

if ~isfield(header,'segments')
    return
end
    
[~,~,EEG_format,EEG_fileS] = lab_filename(cfg.EEG_file);
EEG_fileS = EEG_fileS(1:end-1);

% stitch short segments
disp('Define stitching...')
cfg.STITCH.stitch = {};

lengths = 0;
stitch = {};
Segments = header.segments;
for m = 1:size(Segments,1)
    newlength = (Segments(m,2) - Segments(m,1) + 1) / header.samplingrate;
    stitch{end+1,1} = fullfile(cfg.EEG_filepath,[EEG_fileS num2str(m) '.' EEG_format]); %#ok<AGROW>
    lengths = lengths + newlength;
    if lengths < cfg.STITCH.length | cfg.STITCH.stitchall == true
        stitch{end,2} = 1;
    else
        stitch{end,2} = 0;
        lengths = 0;
    end
end
if stitch{end,2} == 1
    stitch{end,2} = 0;
    tmp = cell2mat(stitch(1:end-1,2));
    tmp2 = find(tmp==0,1,'last');
    if ~isempty(tmp2)
        tmp(tmp2) = 1;
        stitch(1:end-1,2) = num2cell(tmp);
    end
    clearvars tmp tmp2
end
findex = 1;
for j = 1:size(stitch,1)
    fileindex{j,1} = num2str(findex,'%02d'); %#ok<AGROW>
    if stitch{j,2} == 0
        findex = findex+1;
    end
end
stitch = cat(2,stitch,fileindex);
cfg.STITCH.stitch = cat(1,cfg.STITCH.stitch,stitch);
