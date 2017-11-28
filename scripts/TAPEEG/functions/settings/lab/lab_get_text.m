function Text = lab_get_text(Text,Title)

if ~exist('Title','var')
    Title = 'Edit Text';
end
if ~exist('Text','var')
    Text = '';
end
TMP.Text = Text;

Prompt = cell(0,2);
Formats = [];

Prompt(end+1,:) = {Title,Text};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'text';
Formats(end,1).size = 100;

[TMP,Cancelled] = inputsdlg(Prompt,Title,Formats,TMP);
if Cancelled == 1
    Text = '';
    return
else
    Text =TMP.Text;
end