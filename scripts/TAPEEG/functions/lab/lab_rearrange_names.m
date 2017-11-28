function Order = lab_rearrange_names(Names,Title)

if ~exist('Title','var')
    Title = 'Re-arrange Values';
end

Prompt = cell(0,2);
Formats = [];
for i = 1:length(Names)
    Prompt(end+1,:) = {[num2str(i) '. Value'],['Value' num2str(i)]}; %#ok<AGROW>
    Formats(end+1,1).type = 'list'; %#ok<AGROW>
    Formats(end,1).style = 'popupmenu';
    Formats(end,1).items = Names;
    Formats(end,1).callback = {@switch_item,'@ALL','@ALL',i};
    Answer.(['Value' num2str(i)]) = i;
end

[Answer,Cancelled] = inputsdlg(Prompt,Title,Formats,Answer);
if Cancelled == 1
    Order = 1:length(Names);
else
    Order = zeros(1,length(Names));
    for i = 1:length(Names)
        Order(1,i) = Answer.(['Value' num2str(i)]);
    end
end

    function settings = switch_item(settings,value)
        Values = zeros(1,length(Names));
        for j = 1:length(Names)
            Values(1,j) = settings.(['Value' num2str(j)]);
        end
        ActValue = settings.(['Value' num2str(value)]);
        MissValue = setdiff(1:length(Names),Values);
        if ~isempty(MissValue)
            Values(value) = 0;
            tmp = find(Values==ActValue,1,'first');
            if ~isempty(tmp)
                settings.(['Value' num2str(tmp)]) = MissValue;
            end
        end
    end

end