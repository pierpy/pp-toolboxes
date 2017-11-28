function settings = lab_get_DTF(settings,data,header)

if ~exist('data','var')
    data = [];
end
if ~exist('header','var')
    header = [];
end

if exist('settings','var') & isfield(settings,'measure') & ~any(strcmp(settings.measure,'DTF'))
    settings.DTF = [];
    return
end

if ~exist('settings','var') | ~isfield(settings,'DTF') | ~isfield(settings.DTF,'markers')
    settings.DTF.lag = 0;
    if size(data,2) >= 2000 | isempty(data)
        settings.DTF.length = 2000;
    else
        settings.DTF.length = size(data,2);
    end
    settings.DTF.order = 3;
    settings.DTF.sigtest = false;
    settings.DTF.siglevel = 0.05;
    settings.DTF.shufftimes = 1000;
    if ~isfield(settings,'freqs')
        settings.DTF.lowfreq = 8;
        settings.DTF.highfreq = 13;
    end
end

Prompt = cell(0,2);
Formats = {};

Prompt(end+1,:) = {'Lag','lag'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'integer';
Formats(end,1).limits = [0 inf];
Formats(end,1).size = 60;

Prompt(end+1,:) = {'Length','length'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'integer';
Formats(end,1).limits = [0 inf];
Formats(end,1).size = 60;

Formats(end+1,1).type = 'none';
Formats(end,1).span = [1 2];

Prompt(end+1,:) = {'Order','order'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'integer';
Formats(end,1).limits = [0 inf];
Formats(end,1).size = 60;

Prompt(end+1,:) = {'Optimal order','Button'};
Formats(end+1,1).type = 'button';
Formats(end,1).style = 'pushbutton';
Formats(end,1).size = [130 25];
Formats(end,1).callback = {@calc_order,'@ALL','@ALL',header,data};

if isfield(settings.DTF,'lowfreq')
    Formats(end+1,1).type = 'none';
    Formats(end,1).span = [1 2];
    
    Prompt(end+1,:) = {'Low Freq','lowfreq'};
    Formats(end+1,1).type = 'edit';
    Formats(end,1).format = 'float';
    Formats(end,1).limits = [0 inf];
    Formats(end,1).size = 60;
    
    Prompt(end+1,:) = {'High Freq','highfreq'};
    Formats(end+1,1).type = 'edit';
    Formats(end,1).format = 'float';
    Formats(end,1).limits = [0 inf];
    Formats(end,1).size = 60;
end

Formats(end+1,1).type = 'none';
Formats(end,1).span = [1 2];

Prompt(end+1,:) = {'Input significance test','sigtest'};
Formats(end+1,1).type = 'check';
Formats(end,1).span = [1 2];

Prompt(end+1,:) = {'Shuffling Times','shufftimes'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'float';
Formats(end,1).limits = [0 inf];
Formats(end,1).size = 60;

Prompt(end+1,:) = {'Significance Level','siglevel'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'float';
Formats(end,1).limits = [0 inf];
Formats(end,1).size = 60;

Formats(end+1,1).type = 'none';
Formats(end,1).span = [1 2];

[settings.DTF,Cancelled] = inputsdlg(Prompt,'DTF',Formats,settings.DTF);
if isempty(settings.DTF) | Cancelled == 1
    settings.DTF = [];
    return
end

end

function settings = calc_order(settings,header,data)

if isempty(data)
    settings.order = [];
    return
end

startpnt = [];
if isfield(settings,'markers') & ~isempty(settings.markers) & isfield(header,'events') & ~isempty(header.events)
    if strcmp(settings.markers{1},'Start')
        startpnt = cat(1,startpnt,1);
    else
        tmp = strcmp(header.events.TYP,settings.markers{1});
        if ~isempty(tmp)
            startpnt = cat(1,startpnt,header.events.POS(1,tmp)');
        end
        clearvars tmp
    end
end
if isempty(startpnt)
    startpnt = 1;
end
startpnt = startpnt((startpnt+settings.length-1) < size(data,2));
if isempty(startpnt)
    startpnt = 1;
    settings.length = size(data,2);
end
ts =  data(:,startpnt(1):(startpnt(1) + settings.length-1))';

sbc_error = zeros(1,20);
for i = 1:20
    [~,~,~,SBC] = arfit(ts,i,i);
    sbc_error(i) = SBC;
end

settings.order = find(sbc_error == min(sbc_error),1);

end
