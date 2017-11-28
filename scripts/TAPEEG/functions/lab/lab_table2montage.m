% create montage structure from table with montage information
%
% montage = lab_table2montage(table)
%
% written by F. Hatz 2013

function montage = lab_table2montage(table)

if ~isempty(table)
    montage.name = 'Montage';
    montage.label = table(:,1);
    for i = 1:size(table,1)
        if strcmp(table{i,2}(1:1),'A')
            montage.chans{i,1} = str2num(table{i,2}(2:end));
            montage.chans{i,2} = 1;
        else
            montage.chans{i,1} = str2num(table{i,2});
            montage.chans{i,2} = 0;
        end
        if strcmp(table{i,3}(1:1),'A')
            montage.chans{i,3} = str2num(table{i,3}(2:end));
            montage.chans{i,4} = 1;
        else
            montage.chans{i,3} = str2num(table{i,3});
            montage.chans{i,4} = 0;
        end
    end
    tmp = [montage.chans(:,1:2);montage.chans(:,3:4)];
    tmp1 = max(cellfun(@max,tmp(:,1)));
    tmp2 = max(cell2mat(tmp(:,2)) * tmp1);
    montage.numchans = max(tmp1,tmp2);
    clearvars tmp
else
    montage = [];
end