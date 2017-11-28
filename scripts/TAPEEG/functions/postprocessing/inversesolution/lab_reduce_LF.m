function LF = lab_reduce_LF(LF,includelocs,includespi)

if ~exist('LF','var')
   [tmp1,tmp2] = uigetfile('*.bin','Select Leadfield-File');
   if isnumeric(tmp1)
       LF = [];
       return
   end
   LF = fullfile(tmp2,tmp1);
   clearvars tmp1 tmp2
end

if ischar(LF)
    if exist(LF,'file')
        LF_file = LF;
        LF_fileS = LF_file(1:end-7);
        if exist([LF_fileS '.els'],'file')
            locs = lab_read_locs([LF_fileS '.els']);
            if ~isempty(locs)
                numlocs = size(locs.x,2);
            else
                numlocs = [];
            end
        else
            numlocs = [];
        end
        if exist([LF_fileS '.spi'],'file')
            slocs = lab_read_spi([LF_fileS '.spi']);
            if ~isempty(slocs)
                numspi = size(slocs.x,2);
            else
                numspi = [];
            end
        else
            numspi = [];
        end
        LF = lab_read_LFbin(LF_file,numlocs,numspi);
    else
        return
    end
end

if ~exist('includelocs','var') & ~exist('includespi','var')
    if ~exist('locs','var') | isempty(locs) | isempty(locs.labels{1})
        for i = 1:size(LF,1)
            labels{1,i} = ['Channel_' num2str(i)];
        end
    else
        labels = locs.labels;
    end
    if ~exist('slocs','var') | isempty(slocs) | isempty(slocs.labels{1})
        for i = 1:size(LF,3)
            slabels{1,i} = ['SP_' num2str(i)];
        end
    else
        slabels = slocs.labels;
    end
    
    settings.locs = 1:size(LF,1);
    settings.spi = 1:size(LF,3);
    
    Prompt = cell(0,2);
    Formats = [];
    
    Prompt(end+1,:) = {'Electrodes','locs'};
    Formats(end+1,1).type = 'list';
    Formats(end,1).style = 'listbox';
    Formats(end,1).items = labels;
    Formats(end,1).limits = [0 4];
    Formats(end,1).size = [100 300];
    
    Prompt(end+1,:) = {'Solutionpoints','spi'};
    Formats(end+1,1).type = 'list';
    Formats(end,1).style = 'listbox';
    Formats(end,1).items = slabels;
    Formats(end,1).limits = [0 4];
    Formats(end,1).size = [100 300];
    
    [settings,Cancelled] = inputsdlg(Prompt,'Reduce Leadfield',Formats,settings,2);
    if Cancelled == 1
        LF = [];
        return
    else
        includelocs = settings.locs;
        includespi = settings.spi;
    end
end

if ~exist('includelocs','var') | isempty(includelocs)
    includelocs = 1:size(LF,1);
end
if ~exist('includespi','var') | isempty(includespi)
    includespi = 1:size(LF,3);
end
LF = LF(includelocs,:,includespi);

if exist('LF_fileS','var')
    LF_fileS = [LF_fileS '_reduced'];
    lab_write_LFbin([LF_fileS '_LF.bin'],LF);
    if exist('locs','var') & ~isempty(locs)
        exclude = setdiff(1:numlocs,includelocs);
        locs = lab_reduce_locs(locs,exclude);
        lab_write_els([LF_fileS '.els'],locs);
    end
    if exist('slocs','var') & ~isempty(slocs)
        exclude = setdiff(1:numspi,includespi);
        slocs = lab_reduce_locs(slocs,exclude);
        lab_write_spi([LF_fileS '.spi'],slocs);
    end
    LF = [LF_fileS '_LF.bin'];
end
    
end