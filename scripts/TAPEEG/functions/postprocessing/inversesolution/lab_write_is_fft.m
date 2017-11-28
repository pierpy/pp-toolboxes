% Write results of lab_inversesolution_fft
%
% written by F. Hatz 2012

function cfg = lab_write_is_fft(ISresult,cfg)

[~,~,~,cfg.Output_fileS] = lab_filename(cfg.Output_file);

for i = 1:size(cfg.ISFFT.eformat,1)
    cfg.ISspect_file = [cfg.Output_fileS '.' cfg.ISFFT.eformat{i,1}];
    filename = fullfile(cfg.ISspect_filepath,cfg.ISspect_file);
    if strcmp(cfg.ISFFT.eformat{i,1},'edf') & size(ISresult.data,1) > 200
        skipwriting = 1;
        disp(['   skip writing of ' cfg.ISFFT.eformat{i,1} ' (to many channels)']);
    elseif ~strcmp(cfg.ISFFT.eformat{i,1},'ris') & size(ISresult.data,1) > 400
        skipwriting = 1;
        disp(['   skip writing of ' cfg.ISFFT.eformat{i,1} ' (to many channels)']);
    else
        skipwriting = 0;
    end
    if skipwriting == 0
        if strcmp(cfg.ISFFT.eformat{i,1},'ris')
            [~,cfg.ISFFT] = lab_save_data(ISresult,[],cfg.ISFFT.eformat{i,1},filename,cfg.ISFFT,cfg);
        else
            [~,cfg.ISFFT] = lab_save_data(ISresult.data,ISresult,cfg.ISFFT.eformat{i,1},filename,cfg.ISFFT,cfg);
        end
    end
end
tmp = cd(cfg.ISspect_filepath);
if isfield(ISresult,'channels')
    channels = cat(1,num2cell(1:ISresult.numchannels),cellstr(ISresult.channels)');
    channels = cat(1,channels,ISresult.locs.labels);
    lab_write_xls([cfg.Output_fileS '_chans.xls'],channels');
end
if isfield(ISresult,'spectrum') && size(ISresult.data,2) < 257
    xlsresult = [cellstr('Freqs') ISresult.spectrum];
    xlsresult = cat(1,xlsresult,cat(2,cellstr(ISresult.channels),num2cell(ISresult.data)));
    if size(xlsresult,2) > 255
        fileout = [cfg.Output_fileS '.xlsx'];
    else
        fileout = [cfg.Output_fileS '.xls'];
    end
    lab_write_xls(fileout,xlsresult);
end
cd(tmp)
clearvars tmp channels
