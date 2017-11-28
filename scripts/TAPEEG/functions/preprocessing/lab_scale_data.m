% Scale data
%
% [data,header,cfg] = lab_scale_data(data,header,cfg)
%
% data     = matrix (chans x timeframes)
% header   = output of lab_read_data
% cfg      = structure with config (optional)
%
% written by F. Hatz 2012

function [data,header,cfg,skipprocessing] = lab_scale_data(data,header,cfg)
disp('   Scale MEG data')

skipprocessing = 0;

if ~isfield(cfg,'EEG_file') & isfield(header,'EEG_file')
    cfg.EEG_file = header.EEG_file;
    cfg.EEG_filepath = header.EEG_filepath;
end
if ~exist('header','var') | ~isfield(header,'numchannels')
    header.numchannels = size(data,1);
end

if ~isfield(cfg,'SCALE') | ~isfield(cfg.SCALE,'scales')
    [cfg,skipprocessing] = lab_set_scale_data(cfg,header);
end

if skipprocessing == 0
    for i = 1:size(cfg.SCALE.scales,1)
        data(cfg.SCALE.scales(i,1):cfg.SCALE.scales(i,2),:) = data(cfg.SCALE.scales(i,1):cfg.SCALE.scales(i,2),:) * cfg.SCALE.scales(i,3);
    end
    
    % Write verbose file
    if isfield(cfg,'EEG_file')
        Verbose_file=[cfg.EEG_file(1:end-4) '_scale.vrb'];
        fid=fopen(fullfile(cfg.EEG_filepath,Verbose_file),'w');
        fprintf(fid,'Scale data\n');
        fprintf(fid,datestr(now,0));
        fprintf(fid,'\n');
        fprintf(fid,['Input file: ' cfg.EEG_file]);
        fprintf(fid,'\n\n');
        for i = 1:size(cfg.SCALE.scales,1)
            fprintf(fid,['Scale-factor ' num2str(cfg.SCALE.scales(i,1)) '-' num2str(cfg.SCALE.scales(i,2)) ': ' num2str(cfg.SCALE.scales(i,3),'%e')]);
            fprintf(fid,'\n\n');
        end
        fclose(fid);
    end
end

return