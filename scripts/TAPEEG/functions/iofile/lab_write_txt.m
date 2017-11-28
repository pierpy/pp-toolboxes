% Write eeg/meg in ascii (txt)
%
% lab_write_txt(filename,data,header)
%
% written by F. Hatz 2012

function lab_write_txt(filename,data,header)
[~,filepath,~,filename] = lab_filename(filename);
filename = fullfile(filepath,[filename '.txt']);

fid=fopen(filename,'w');
if isnumeric(data)
    if exist('header','var') & isfield(header,'channels')
        for i = 1:size(data,1)
            fprintf(fid,header.channels(i,:));
            fprintf(fid,'\t');
        end
        fprintf(fid,native2unicode([13 10]));
    end
    for i = 1:size(data,2)
        fprintf(fid,'%10.5f\t',data(:,i));
        fprintf(fid,native2unicode([13 10]));
    end
elseif iscell(data)
    for i = 1:size(data,1)
        for j = 1:size(data,2)
            if ischar(data{i,j})
                fprintf(fid,'%s',data{i,j});
            elseif isnumeric(data{i,j})
                fprintf(fid,'%10.5f',data{i,j});
            else
                fprintf(fid,'%s','NaN');
            end
            if j < size(data,2)
                fprintf(fid,'%s\t','');
            end
        end
        fprintf(fid,native2unicode([13 10]));  
    end
end
fclose(fid);

return
