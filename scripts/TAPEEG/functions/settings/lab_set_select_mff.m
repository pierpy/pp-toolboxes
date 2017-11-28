function [settings,skipprocessing] = lab_set_select_mff(settings,header)

skipprocessing = 0;

if ~exist('header','var')
    header = [];
end
if ~exist('settings','var')
    settings = [];
end
if ~isfield(settings,'SELECTMFF') | isempty(settings.SELECTMFF)
    settings.SELECTMFF.number = [];
    settings.SELECTMFF.markerinclude = '';
    settings.SELECTMFF.markerexclude = '';
    settings.SELECTMFF.minlength = [];
    settings.SELECTMFF.minlength = [];
    settings.SELECTMFF.longestsegment = false;
end

if isfield(header,'events') & isfield(header.events,'TYP') & ~isempty(header.events.TYP)
    markerlist = unique(header.events.TYP);
    markerlist = markerlist(:)';
else
    markerlist = [];
end

Prompt = cell(0,3);
Formats = [];

Prompt(end+1,:) = {'Select by number','number',''};
Formats(end+1,1).type = 'check';
Formats(end,1).callback = {@set_number,'@ALL','@ALL',header};

Formats(end+1,1).type = 'none';

Prompt(end+1,:) = {'Include segments with marker(s):','',''};
Formats(end+1,1).type = 'text';

Prompt(end+1,:) = {'','markerinclude',''};
if ~isempty(markerlist)
    Formats(end+1,1).type = 'edit';
    Formats(end,1).format = 'vector';
    Formats(end,1).items = markerlist;
    Formats(end,1).size = 250;
else
    Formats(end+1,1).type = 'edit';
    Formats(end,1).format = 'text';
    Formats(end,1).size = 250;
end

Formats(end+1,1).type = 'none';

Prompt(end+1,:) = {'Exclude segments with marker(s):','',''};
Formats(end+1,1).type = 'text';

Prompt(end+1,:) = {'','markerexclude',''};
if ~isempty(markerlist)
    Formats(end+1,1).type = 'edit';
    Formats(end,1).format = 'vector';
    Formats(end,1).items = markerlist;
    Formats(end,1).size = 250;
else
    Formats(end+1,1).type = 'edit';
    Formats(end,1).format = 'text';
    Formats(end,1).size = 250;
end

Formats(end+1,1).type = 'none';

Prompt(end+1,:) = {'Minimal length','minlength','seconds'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'integer';
Formats(end,1).limits = [0 999999999];
Formats(end,1).size = 40;

Prompt(end+1,:) = {'Maximal length','maxlength','seconds'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'integer';
Formats(end,1).limits = [0 999999999];
Formats(end,1).size = 40;

Formats(end+1,1).type = 'none';

Prompt(end+1,:) = {'Take longest segment','longestsegment',''};
Formats(end+1,1).type = 'check';
Formats(end,1).format = 'input';
Formats(end,1).callback = {@take_longest,'@ALL','@ALL'};


[settings.SELECTMFF,Cancelled] = inputsdlg(Prompt,'Segments selection',Formats,settings.SELECTMFF);
if Cancelled == 1
    settings.SELECTMFF = [];
    return
end

end

function settings = set_number(settings,header)
    if isfield(header,'events') & isfield(header.events,'POS') & ~isempty(header.events.POS)
        markerstart = header.events.POS(1,ismember(header.events.TYP,'SegStart')==1);
        markerend = header.events.POS(1,ismember(header.events.TYP,'SegStop')==1);
        if size(markerstart,2) == size(markerend,2) & size(markerstart,2) > 0
            header.segments = [markerstart' markerend'];
            for i = 1:size(header.segments,1)
                seglength = (header.segments(i,2) - header.segments(i,1) + 1) / header.samplingrate;
                seglist{i} = ['Segment' num2str(i) ' (' num2str(seglength) ' sec)']; %#ok<AGROW>
            end
        else
            settings.number = [];
            disp('   No valid segments found')
            return
        end
    else
        for i = 1:20
            seglist{i} = ['Segment' num2str(i)]; %#ok<AGROW>
        end
    end
    Prompt = {'Segments','number'};
    Formats.type = 'list';
    Formats(end,1).style = 'listbox';
    Formats(end,1).limits = [0 4];
    Formats.size = [190 170];
    Formats.items = seglist;
    [settings,Cancelled] = inputsdlg(Prompt,'Select Numbers',Formats,settings);
    if Cancelled == 1
        settings.number = [];
    end
end

function settings = take_longest(settings)
    if settings.longestsegment == false
        settings.number = [];
        settings.markerinclude = '';
        settings.markerexclude = '';
        settings.minlength = [];
        settings.minlength = [];
        settings.longestsegment = true;
    else
        settings.longestsegment = false;
    end
end