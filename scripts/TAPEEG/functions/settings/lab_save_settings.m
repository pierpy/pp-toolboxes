function lab_save_settings(cfg,calc)

if (isfield(cfg,'MAIN') & isfield(cfg.MAIN,'auto') & cfg.MAIN.auto == 1) | ...
        (isfield(cfg,'editsettings') & cfg.editsettings == false)
    fid = fopen(fullfile(cfg.settings_path,'settings.vrb'),'a');
    fprintf(fid,'Settings file\n');
    fprintf(fid,datestr(now,0));
    fprintf(fid,'\n');
    fprintf(fid,strrep(cfg.settings_path,'\','/'));
    fprintf(fid,'/');
    fprintf(fid,cfg.settings_file);
    fprintf(fid,'\n\n');
    fclose(fid);
    return
end

if ~isfield(cfg,'settings_file')
    cfg.settings_file = 'settings.mat';
end
if ~isfield(cfg,'settings_path')
    cfg.settings_path = pwd;
end

filename = fullfile(cfg.settings_path,cfg.settings_file);

if isfield(cfg,'settings_file')
    cfg = rmfield(cfg,'settings_file');
end
if isfield(cfg,'settings_path')
    cfg = rmfield(cfg,'settings_path');
end
if isfield(cfg,'EEG_file')
    cfg = rmfield(cfg,'EEG_file');
end
if isfield(cfg,'EEG_filepath')
    cfg = rmfield(cfg,'EEG_filepath');
end
if isfield(cfg,'MAIN') & isfield(cfg.MAIN,'filelist') & exist('calc','var')
    cfg.MAIN.filelist = union(cfg.MAIN.filelist,calc.Filelist_done);
    cfg.MAIN.filelist = cfg.MAIN.filelist(:)';
elseif exist('calc','var')
    cfg.MAIN.filelist = calc.Filelist_done;
else
    cfg.MAIN.filelist = [];
end
if isfield(cfg,'listold')
    cfg = rmfield(cfg,'listold');
end
if isfield(cfg,'MAIN') & isfield(cfg.MAIN,'log_path')
    cfg.MAIN = rmfield(cfg.MAIN,'log_path');
end
if isfield(cfg,'STITCH') & isfield(cfg.STITCH,'stitch')
    cfg.STITCH = rmfield(cfg.STITCH,'stitch');
end
if isfield(cfg,'lastsegment')
    cfg = rmfield(cfg,'lastsegment');
end

cfg.timestamp = date;

if exist(filename,'file')
    delete(filename)
end
save(filename,'cfg','-v7.3');