function Marker = lab_get_marker(Marker)

if ischar(Marker) & strcmp(Marker,'**edit**')
    TMP.Marker = '';
    Prompt = cell(0,2);
    Formats = [];
    Prompt(end+1,:) = {'Marker','Marker'};
    Formats(end+1,1).type = 'edit';
    Formats(end,1).format = 'text';
    Formats(end,1).size = 100;
    [TMP,Cancelled] = inputsdlg(Prompt,'Marker',Formats,TMP);
    if Cancelled == 1
        Marker = '';
        return
    else
        Marker =TMP.Marker;
    end
elseif iscell(Marker) & any(strcmp(Marker,'**edit**'))
    TMP.Marker = Marker(~strcmp(Marker,'**edit**'));
    Prompt = cell(0,2);
    Formats = [];
    Prompt(end+1,:) = {'Marker','Marker'};
    Formats(end+1,1).type = 'edit';
    Formats(end,1).format = 'vector';
    Formats(end,1).limits = [0 1];
    Formats(end,1).size = 150;
    [TMP,Cancelled] = inputsdlg(Prompt,'Marker',Formats,TMP);
    if Cancelled == 1
        Marker = {};
        return
    else
        Marker =TMP.Marker;
    end
end