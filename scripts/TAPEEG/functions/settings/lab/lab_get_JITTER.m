function settings = lab_get_JITTER(settings,markerOffset)

if ~exist('settings','var') | isempty(settings)
    if markerOffset > 0
        editpre = false;
        settings.prerange = markerOffset;
    else
        editpre = true;
        settings.prerange = 50;
    end
    settings.percentmax = 70;
    settings.range = 80;
    settings.rangej = 20;
    settings.SD = 1.5;
elseif markerOffset > 0
    editpre = false;
    settings.prerange = markerOffset;
else
    editpre = true;
end

Prompt = cell(0,2);
Formats = [];

Prompt(end+1,:) = {'Test for valid channels:',''};
Formats(end+1,1).type = 'text';

if editpre == true
    Prompt(end+1,:) = {'Pre-Marker range','prerange'};
    Formats(end+1,1).type = 'edit';
    Formats(end,1).format = 'integer';
    Formats(end,1).limits = [0 999];
    Formats(end,1).size = 40;
else
    Prompt(end+1,:) = {['Pre-Marker range: ' num2str(settings.prerange)],''};
    Formats(end+1,1).type = 'text';
end

Prompt(end+1,:) = {'Percent of max amplitude','percentmax'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'float';
Formats(end,1).limits = [1 100];
Formats(end,1).size = 30;

Prompt(end+1,:) = {'  ',''};
Formats(end+1,1).type = 'text';

Prompt(end+1,:) = {'Jitter calculation:',''};
Formats(end+1,1).type = 'text';

Prompt(end+1,:) = {'Timeframes to correlate','range'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'integer';
Formats(end,1).limits = [0 9999];
Formats(end,1).size = 40;

Prompt(end+1,:) = {'Maximal shift','rangej'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'integer';
Formats(end,1).limits = [0 9999];
Formats(end,1).size = 40;

Prompt(end+1,:) = {'Standard deviations','SD'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'float';
Formats(end,1).limits = [0 99];
Formats(end,1).size = 30;

[settings,Cancelled] = inputsdlg(Prompt,'Jitter settings',Formats,settings);
if Cancelled == 1
    settings = [];
else
    pause(0.2);
    settings.editpre = editpre;
end

end