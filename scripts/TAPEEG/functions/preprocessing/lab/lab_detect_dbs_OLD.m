function [data,header,cfg,skipprocessing] = lab_detect_dbs(data,header,cfg)

disp('   detect stimulator artifacts')

skipprocessing = 0;

if ~exist('cfg','var') | ~isfield(cfg,'STIM') | ~isfield(cfg.STIM,'DETECT') |...
         ~isfield(cfg.STIM.DETECT,'stimstd') | isempty(cfg.STIM.DETECT.stimstd)
    cfg.STIM.DETECT.stimstd = 2;
end
if ~isfield(cfg.STIM.DETECT,'stimpercent') | isempty(cfg.STIM.DETECT.stimpercent)
    cfg.STIM.DETECT.stimpercent = 30;
end
if ~exist('header','var')
    header= lab_create_header(data);
end
if ~isfield(header,'numdatachannels')
    header.numdatachannels = size(data,1);
end

MARK = false(header.numdatachannels,size(data,2));
for i = 1:header.numdatachannels
    Idx = find(diff(sign(diff(data(i,:))))== -2) + 1;
    Val = data(i,Idx);
    Idx = Idx(Val > (cfg.STIM.DETECT.stimstd * std(Val)));
    MARK(i,Idx) = true;
end
MARK = sum(MARK,1);
MARK = find(MARK > cfg.STIM.DETECT.stimpercent*header.numdatachannels/100);
DIFF = median(diff(MARK));
DIFF = floor(DIFF/2);
if length(MARK) < 3
    disp('    no high-voltage stimulator found')
    skipprocessing = 1;
    return
end
if MARK(1) < DIFF+1
    MARK = MARK(2:end);
end
if MARK(end) + DIFF > size(data,2)
    MARK = MARK(1:end-1);
end
MARK = MARK(find(diff(MARK) >= 2*DIFF-1)+1);

MARK2 = [];
for i = -DIFF:DIFF
    MARK2 = cat(1,MARK2,MARK+i);
end
stim = data(:,MARK2(:)');
stim = reshape(stim,[size(stim,1) size(MARK2,1),size(MARK2,2)]);
stim = mean(stim,3);
stim2 = stim;
for i = 1:header.numdatachannels
    tmp = stim(i,1):(stim(i,end)-stim(i,1))/(size(stim,2)-1):stim(i,end);
    if ~isempty(tmp)
        stim2(i,:) = stim2(i,:) - tmp;
    end
end
clearvars tmp
header.STIM.raw = stim;
header.STIM.corr = stim2;
header.STIM.marker = MARK;
events.POS = int64(MARK(:)');
events.DUR = int64(ones(1,length(MARK)));
events.OFF = int64(zeros(1,length(MARK)));
events.TYP = repmat({'STIM'},1,length(MARK));
header = lab_mix_markers(header,events);

end

    
    