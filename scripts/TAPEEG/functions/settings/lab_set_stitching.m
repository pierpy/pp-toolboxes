function [cfg,skipprocessing] = lab_set_stitching(cfg)

skipprocessing = 0;

if ~exist('cfg','var') | ~isfield(cfg,'STITCH') | ~isfield(cfg.STITCH,'length')
    cfg.STITCH.length = 35;
    cfg.STITCH.window = 1;
    cfg.STITCH.stitchall = false;
end

Prompt = cell(0,2);
Formats = {};

Prompt(end+1,:) = {'Stitch segments < ...seconds (0=off)','length'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'float';
Formats(end,1).limits = [0 inf];
Formats(end,1).callback = {@do_length,'@ALL','@ALL'};
Formats(end,1).size = 30;

Prompt(end+1,:) = {'Stitch folder' 'stitchall'};
Formats(end+1,1).type = 'check';
Formats(end,1).format = 'input';
Formats(end,1).callback = {@do_stitchall,'@ALL','@ALL'};

Prompt(end+1,:) = {'overlap window (seconds)','window'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'float';
Formats(end,1).limits = [0 inf];
Formats(end,1).size = 30;

[cfg.STITCH,Cancelled] = inputsdlg(Prompt,'Stitch setting',Formats,cfg.STITCH,2);
if Cancelled == 1
    cfg.STITCH = [];
    skipprocessing = 1;
    return
end

    function settings = do_stitchall(settings)
        if settings.stitchall == false
            settings.stitchall = true;
            settings.length = [];
        else
            settings.stitchall = false;
        end
    end
    
    function settings = do_length(settings)
        if settings.stitchall == true & ~isempty(settings.length)
            settings.stitchall = false;
        end
    end
 
end