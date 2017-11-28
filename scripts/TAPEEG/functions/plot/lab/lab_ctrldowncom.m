% Helper file for lab_plot_eeg
%
% Written by F. Hatz 2013

function lab_ctrldowncom(quick_evtmk,quick_evtrm)

ax1 = findobj('tag','backeeg','parent',gcbf);
tmppos = get(ax1, 'currentpoint');


% Store "UserData" and "temppos" variables to "g" structures.             
g=get(gcbf, 'UserData');
g.tmppos=tmppos;

g.quick_evtmk=quick_evtmk;
g.quick_evtrm=quick_evtrm;

if isfield(g,'events')
    g.nevents = size(g.events,2);
else
    g.nevents = 0;
end

g.datasize = g.frames;
numchannels = size(g.datastd,1);

if ~isfield(g.events, 'index');
    for i=1:length(g.events);
        g.events(i).index=i;
        g.events(i).proc='none';
    end
end

g.eventedit.WinStartPnt=g.time*g.srate;
g.eventedit.EpochIndex=1;
g.eventedit.PosLat=round(g.tmppos(1,1)+g.eventedit.WinStartPnt);


% Identify selected channel.
% By default use the 'Eelec' tex display of eegplot.
tmpChanIndex=strmatch(get(findobj(gcf,'Tag','Eelec'),'string'),{g.eloc_file.labels},'exact');
if length(tmpChanIndex)==1;
    g.eventedit.ChanIndex=tmpChanIndex;
else
    % Otherwise calculate ChanIndex from tmppos.
    nWin=(g.chans-g.dispchans)+1;
    stepWin=1/nWin;
    if g.dispchans==g.chans;
        curWinStrt=0;
    else
        curWinStrt=floor((1-get(findobj('tag','eegslider'),'value'))/stepWin);
    end
    curWinEnd=curWinStrt+g.dispchans;

    YIndex=floor((tmppos(1,2)/(1/(g.dispchans+1)))-.5);
    g.eventedit.ChanIndex=(curWinEnd-YIndex);

    if g.eventedit.ChanIndex==0;
        g.eventedit.ChanIndex=1;
    end
    if g.eventedit.ChanIndex>numchannels;
        g.eventedit.ChanIndex=numchannels;
    end
end
clear tmpChanIndex

% Check for event selection (if events already exist in dataset).
if isfield(g,'events') & ~isempty(g.events);

    % Check for event selection (within +/-20 points of button press).
    if isfield(g.eventedit, 'SelEventStruct');
        g.eventedit=rmfield(g.eventedit,'SelEventStruct');
    end
    j=0;
    for i=1:length(g.events);
        if ~isfield(g.events,'duration') | isempty(g.events(i).duration)
            g.events(i).duration = 1;
        end
        if g.eventedit.PosLat>(g.events(i).latency-20) & g.eventedit.PosLat<(g.events(i).latency+g.events(i).duration+19);
            j=j+1;
            g.eventedit.SelEventStruct(j).index=i;
            g.eventedit.SelEventStruct(j).dist=abs(g.events(i).latency-g.eventedit.PosLat);
            g.eventedit.SelEventStruct(j).type=g.events(i).type;
            g.eventedit.SelEventStruct(j).latency=g.events(i).latency;
            g.eventedit.SelEventStruct(j).duration=g.events(i).duration;
        end
    end
end

if g.eventedit.PosLat < g.eventedit.WinStartPnt
    if isfield(g,'eloc_file') & isfield(g.eloc_file,'badchan')
        if g.eloc_file(1,g.eventedit.ChanIndex).badchan == 0
            g.eloc_file(1,g.eventedit.ChanIndex).badchan = 1;
        else
            g.eloc_file(1,g.eventedit.ChanIndex).badchan = 0;
        end
    end
elseif isfield(g.eventedit,'SelEventStruct')
    tmp = g.eventedit.SelEventStruct;
    index = g.eventedit.SelEventStruct.index;
    for i = 1:length(tmp)
        events(i).Name = tmp(i).type;
        events(i).Position = tmp(i).latency+1;
        events(i).Duration = tmp(i).duration;
        events(i).Delete = false;
    end
    [events,skipprocessing] = inputsdlg(events,'Edit event');
    if skipprocessing == 1
        return
    end
    deleteflag = [];
    for i = 1:length(index)
        if events(i).Delete == false
            g.events(1,index(i)).type = events(i).Name;
            g.events(1,index(i)).latency = events(i).Position-1;
            g.events(1,index(i)).duration = events(i).Duration;
            if ~isfield(g, 'eventupdate');
                updateindex=1;
            else
                updateindex=length(g.eventupdate)+1;
            end
            g.eventupdate(updateindex).latency = g.events(index(i)).latency;
            g.eventupdate(updateindex).type = g.events(index(i)).type;
            g.eventupdate(updateindex).proc = 'edit';
            g.eventupdate(updateindex).index = g.events(index(i)).index;
        else
            deleteflag = [deleteflag index(i)];
        end
    end
    if ~isempty(deleteflag)
        for i = 1:length(deleteflag)
            if ~isfield(g, 'eventupdate');
                updateindex=1;
            else
                updateindex=length(g.eventupdate)+1;
            end
            g.eventupdate(updateindex).latency = [];
            g.eventupdate(updateindex).type = [];
            g.eventupdate(updateindex).proc = 'clear';
            g.eventupdate(updateindex).index = g.events(deleteflag(i)).index;
        end
        g.events(deleteflag)=[];
    end
else
    events.Name = 'New';
    events.Position = g.eventedit.PosLat+1;
    events.Duration = 1;
    [events,skipprocessing] = inputsdlg(events,'Edit event');
    if skipprocessing == 1
        return
    end
    if ~isfield(g, 'newindex');
        g.newindex=g.nevents+1;
    else
        g.newindex=g.newindex+1;
    end
    g.events(1,end+1).type = events.Name;
    g.events(1,end).latency = events.Position-1;
    g.events(1,end).duration = events.Duration;
    g.events(1,end).proc='new';
    g.events(1,end).index = length(g.events);
    for i = 1:length(g.events)
        tmp(i) = g.events(1,i).latency;
    end
    [~,tmp] = sort(tmp);
    g.events = g.events(tmp);
    
    if ~isfield(g, 'eventupdate');
        updateindex=1;
    else
        updateindex=length(g.eventupdate)+1;
    end
    g.eventupdate(updateindex).latency=events.Duration;
    g.eventupdate(updateindex).type=events.Name;
    g.eventupdate(updateindex).proc='new';
    g.eventupdate(updateindex).index=g.newindex;
end

for i = 1:length(g.events)
    EventType = g.events(i).type;
    if isfield(g, 'eventtypes');
        if ~any(strcmp(EventType, g.eventtypes));
            eventtypesN=length(g.eventtypes)+1;
            g.eventtypes{eventtypesN} = EventType;
            tmp = lines;
            tmp2 = mod(eventtypesN,64);
            g.eventtypecolors{eventtypesN} = tmp(tmp2,:);
            clearvars tmp tmp2
            g.eventtypestyle{eventtypesN} = '-';
            g.eventtypewidths(eventtypesN) = 1;
        end
    else
        eventtypesN=1;
        g.eventtypes{eventtypesN} = EventType;
        g.eventtypecolors{eventtypesN} = 'k';
        g.eventtypestyle{eventtypesN} = '-';
        g.eventtypewidths(eventtypesN) = 1;
        g.plotevent='on';
    end
end

if isfield(g, 'eventcolors');
    fields={'eventcolors', 'eventstyle', 'eventwidths', 'eventlatencies', 'eventlatencyend'};
    g=rmfield(g,fields);
end

if isempty(g.events);
    g.eventcolors=[];
    g.eventstyle=[];
    g.eventwidths=[];
    g.eventlatencies=[];
    g.eventlatencyend=[];    
else
    for i=1:length(g.events);
        eventtypeindex=find(strcmp(g.eventtypes,g.events(i).type));
        g.eventcolors{i}=g.eventtypecolors{eventtypeindex};
        g.eventstyle{i}=g.eventtypestyle{eventtypeindex};
        g.eventwidths(i)=g.eventtypewidths(eventtypeindex);
        g.eventlatencies(i)=g.events(i).latency;
        g.eventlatencyend(i)=g.events(i).latency+g.eventwidths(i)+g.events(i).duration;    
    end
end

g = rmfield(g, 'eventedit');
set(findobj('tag', 'EEGPLOT'), 'UserData', g);
assignin('base','EVENTS',g.events);
assignin('base','BAD',g.eloc_file);
%eegplot('defmotioncom',gcbf);
eegplot('drawp', 0);

% Call event edit UI.


