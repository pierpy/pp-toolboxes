function Number = lab_get_number(Number,Title)

if ~exist('Title','var')
    Title = 'Edit Number';
end
if ~exist('Number','var')
    Number = [];
end
TMP.Number = Number;

Prompt = cell(0,2);
Formats = [];

Prompt(end+1,:) = {Title,'Number'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'integer';
Formats(end,1).limits = [0 9999];
Formats(end,1).size = 30;

[TMP,Cancelled] = inputsdlg(Prompt,Title,Formats,TMP);
if Cancelled == 1
    Number = [];
else
    Number =TMP.Number;
end