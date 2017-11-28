function [Parts,cfg] = lab_define_filereading(header,cfg)
    
disp('Define File reading')

if ~exist('header','var') | ~isfield(header,'numtimeframes') | isempty(header.numtimeframes)
    Parts = [];
    disp('   Defining file reading not possible, no header information')
    return
else
    Parts = [1 header.numtimeframes];
end
if ~isfield(header.events,'TYP')
    domff = false;
else
    markerstart = header.events.POS(1,ismember(header.events.TYP,'SegStart')==1);
    markerend = header.events.POS(1,ismember(header.events.TYP,'SegStop')==1);
    if isempty(markerstart) | length(markerstart) ~= length(markerend)
        domff = false;
    else
        Parts = [markerstart' markerend'];
        domff = true;
    end
end
if ~exist('cfg','var')
    cfg = [];
end
if ~isfield(cfg,'FILEREADING')
    [cfg,skipprocessing] = lab_set_define_filereading(cfg,header);
    pause(0.2);
    if skipprocessing == 1
        Parts = [];
        return
    end
end

if domff == true & isfield(cfg.FILEREADING,'SELECTMFF') & ~isempty(cfg.FILEREADING.SELECTMFF)
    % define mff reading
    flag = true(1,size(Parts,1));
    
    if isfield(cfg.FILEREADING.SELECTMFF,'number') & ~isempty(cfg.FILEREADING.SELECTMFF.number)
        disp('   Select Parts by number')
        flag = false(1,size(Parts,1));
        flag(cfg.FILEREADING.SELECTMFF.number) = true;
    end
    
    if isfield(cfg.FILEREADING.SELECTMFF,'markerinclude') & ~isempty(cfg.FILEREADING.SELECTMFF.markerinclude)
        disp('   Select Parts by markers to include')
        for i = 1:size(Parts,1)
            idx = intersect(find(header.events.POS >= int64(Parts(i,1))),find(header.events.POS <= int64(Parts(i,2))));
            markers = header.events.TYP(idx);
            if isempty(intersect(markers,cfg.FILEREADING.SELECTMFF.markerinclude))
                flag(i) = false;
            end
        end
    end
    
    if isfield(cfg.FILEREADING.SELECTMFF,'markerexclude') & ~isempty(cfg.FILEREADING.SELECTMFF.markerexclude)
        disp('   Select Parts by markers to exclude')
        for i = 1:size(Parts,1)
            idx = intersect(find(header.events.POS >= int64(Parts(i,1))),find(header.events.POS <= int64(Parts(i,2))));
            markers = header.events.TYP(idx);
            if ~isempty(intersect(markers,cfg.FILEREADING.SELECTMFF.markerexclude))
                flag(i) = false;
            end
        end
    end
    
    if isfield(cfg.FILEREADING.SELECTMFF,'minlength') & ~isempty(cfg.FILEREADING.SELECTMFF.minlength)
        disp('   Select Parts by minimal length')
        for i = 1:size(Parts,1)
            if (Parts(i,2) - Parts(i,1) + 1) < cfg.FILEREADING.SELECTMFF.minlength*header.samplingrate
                flag(i) = false;
            end
        end
    end
    
    if isfield(cfg.FILEREADING.SELECTMFF,'maxlength') & ~isempty(cfg.FILEREADING.SELECTMFF.maxlength)
        disp('   Select Parts by maximal length')
        for i = 1:size(Parts,1)
            if (Parts(i,2) - Parts(i,1) + 1) > cfg.FILEREADING.SELECTMFF.maxlength*header.samplingrate
                flag(i) = false;
            end
        end
    end
    
    if max(flag) == 0
        disp('   No Parts match criteria, skip file')
        Parts = [];
        return
    else
        Parts = Parts(flag,:);
    end
    
    if isfield(cfg.FILEREADING.SELECTMFF,'longestsegment') & cfg.FILEREADING.SELECTMFF.longestsegment == true
        disp('   Select longest segment')
        Parts = Parts((Parts(:,2) - Parts(:,1)) == max((Parts(:,2) - Parts(:,1))),:);
    end
end

% define sectional reading
if isfield(cfg.FILEREADING,'sectionalread') & cfg.FILEREADING.sectionalread > 0 & ...
        isfield(header,'sectionalreading') & header.sectionalreading == true
    Parts2 = [];
    for i = 1:size(Parts,1)
        for j = 1:ceil((Parts(i,2)-Parts(i,1)+1) / cfg.FILEREADING.sectionalread)
            Parts2(end+1,1) = Parts(i,1) + (i-1)*cfg.FILEREADING.sectionalread; %#ok<AGROW>
            Parts2(end,2) = Parts(i,1) + i*cfg.FILEREADING.sectionalread;
            if Parts2(end,2) > Parts(i,2)
                Parts2(end,2) = Parts(i,2);
            end
        end
    end
    Parts = Parts2;
    clearvars Parts2 i j
end

end