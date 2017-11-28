function Scaletxt = lab_get_format(Format,Scaletxt)

if ~exist('Scaletxt','var') | isempty(Scaletxt)
    settings.Scaletxt = '1';
else
    settings.Scaletxt = Scaletxt;
end

if (ischar(Format) & strcmp(Format,'txt')) | (iscell(Format) & any(strcmp(Format,'txt')))
    Prompt = cell(0,2);
    Formats = [];
    Prompt(end+1,:) = {'Scaletxt data for txt-export (auto/number/1=off)','Scaletxt'};
    Formats(end+1,1).type = 'edit';
    Formats(end,1).format = 'text';
    Formats(end,1).size = 90; % automatically assign the height
    [settings,Cancelled] = inputsdlg(Prompt,'Text scaling',Formats,settings);
    if Cancelled == 1
        Scaletxt = '1';
    else
        Scaletxt = settings.Scaletxt;
    end
else
    Scaletxt = [];
end