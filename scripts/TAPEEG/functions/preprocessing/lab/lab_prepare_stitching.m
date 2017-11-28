% Helper script for lab_stitching
%
% written by F. Hatz 2012

function lab_prepare_stitching(cfg,calc)

global stitchinfo Segments
    
% stitch short segments
disp('Define stitching...')
for i = 1:size(calc.Filelist,2)
    [EEGfiles{i},EEGpath{i}] = lab_filename(calc.Filelist{1,i}); %#ok<AGROW>
end
tmp = find(strcmp(EEGpath,cfg.EEG_filepath));
EEGpath = EEGpath(tmp);
EEGfiles = EEGfiles(tmp);
stitch = {};
lengths = 0;
for i = 1:length(EEGfiles)
    [~,~,EEG_format,EEG_fileS] = lab_filename(EEGfiles{i});
    SegmentsTmp = {};
    if strcmp(EEG_fileS,cfg.EEG_fileS) | strcmp([EEG_fileS '_S1'],cfg.EEG_fileS)
        SegmentsTmp = Segments(2,:);
    else
        disp(['    read file ' fullfile(EEGpath{i},EEGfiles{i})])
        datatmp = [];
        headertmp = [];
        try %#ok<TRYNC>
            if isfield(cfg,'SEG') & isfield(cfg.SEG,'select') & strcmp(cfg.SEG.select,'Automated')
                [datatmp,headertmp,cfg] = lab_read_data(fullfile(EEGpath{i},EEGfiles{i}),cfg);
            else
                [~,headertmp,cfg] = lab_read_data(fullfile(EEGpath{i},EEGfiles{i}),cfg,true);
            end
        end
        if ~isempty(headertmp)
            if i == length(EEGfiles)
                cfg.lastfilefolder = true;
            else
                cfg.lastfilefolder = false;
            end
            if i == 1
                cfg.firstfilefolder = true;
            else
                cfg.firstfilefolder = false;
            end
            [SegmentsTmp,cfg] = lab_define_segments(datatmp,headertmp,cfg);
            SegmentsTmp = SegmentsTmp(2,:);
            [~,~,EEG_format,EEG_fileS] = lab_filename(cfg.EEG_file);
        end
    end
    for j = 1:size(SegmentsTmp,2)
        if size(SegmentsTmp,2) > 1
            EEG_file = [EEG_fileS '_S' num2str(j) '.' EEG_format];
        else
            EEG_file = [EEG_fileS '.' EEG_format];
        end
        stitch{end+1,1} = fullfile(cfg.EEG_filepath,EEG_file); %#ok<AGROW>
        newlength = SegmentsTmp{1,j}.numtimeframes/SegmentsTmp{1,j}.samplingrate;
        lengths = lengths + newlength;
        if lengths < cfg.STITCH.length | cfg.STITCH.stitchall == true
            stitch{end,2} = 1;
        else
            stitch{end,2} = 0;
            lengths = 0;
        end
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
stitchinfo = cat(1,stitchinfo,stitch);

return