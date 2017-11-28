function cfg = lab_convert_els(cfg,numdatachannels,selected)
fid=fopen(fullfile(cfg.settings_path,'electrodes.els'),'r');
if fid > 0
    head.version = fgets(fid);
    head.number = fgets(fid);
    head.add1 = fgets(fid);
    head.add2 = fgets(fid);
    head.add3 = fgets(fid);
    head.add4 = fgets(fid);
    NumHead = str2num(head.number); %#ok<ST2NM>
    for i = 1:NumHead
        position{i,1} = fgets(fid); %#ok<AGROW>
    end
    fclose(fid);
    if NumHead == numdatachannels
        position = position(selected,:);
        head.number = [num2str(size(position,1)) native2unicode([13 10])];
        if ~isfield(cfg,'ELS_filepath')
            cfg.ELS_filepath = cfg.EEG_filepath;
        end
        if ~isfield(cfg,'ELS_file')
            cfg.ELS_file = [cfg.EEG_file(1:end-4) '.els'];
        else
            cfg.ELS_file = [cfg.ELS_file(1:end-4) '.els'];
        end
        fidout=fopen(fullfile(cfg.ELS_filepath,cfg.ELS_file),'w');
        fprintf(fidout,head.version);
        fprintf(fidout,head.number);
        fprintf(fidout,head.add1);
        fprintf(fidout,head.add2);
        fprintf(fidout,head.number);
        fprintf(fidout,head.add4);
        for i = 1:size(position,1)
            fprintf(fidout,position{i,1});
        end
        fclose(fidout);
    else
        cfg.ELS_file = [];
        warning('electrodes.els not matching')
    end
end
return