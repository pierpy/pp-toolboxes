% Convert montage-structure to text-information (for display / output)

function text = lab_montage2txt(montage,header)

text = sprintf('Name\t\tActive\t\tReference\n');
for j = 1:length(montage)
    text = sprintf([text montage(j).name '\n']);
    for i = 1:size(montage(j).chans,1);
        text = sprintf([text montage(j).label{i,1} '\t\t']);
        active = montage.chans{i,1}(1,1) + (montage.chans{i,2} * header.numdatachannels);
        text = sprintf([text header.channels(active,:) '\t\t']);
        if isnumeric(montage.chans{i,3}) | isempty(montage.chans{i,3})
            reference = montage.chans{i,3} + (montage.chans{i,4} * header.numdatachannels);
            text = sprintf([text num2str(reference)]);
        elseif strcmp(montage.chans{i,3},'AVG')
            text = sprintf([text 'AVG\n']);
        elseif strcmp(montage.chans{i,3},'LAPL')
            text = sprintf([text 'LAPL\n']);
        end
        text = sprintf([text '\n']);
    end
end