% Write electrodes localications in .xyz format
%
% lab_write_xyz(locs,filename)
%
% written by F. Hatz 2012

function lab_write_xyz(filename,locs)
[~,filepath,~,filename] = lab_filename(filename);
filename = fullfile(filepath,[filename '.xyz']);

if isstruct(locs)
    if isfield(locs,'x')
        for i = 1:size(locs.x,2)
            positions{i,1} = num2str(locs.x(1,i),6);
            positions{i,2} = num2str(locs.y(1,i),6);
            positions{i,3} = num2str(locs.z(1,i),6);
            positions{i,4} = locs.labels{1,i};
        end
    elseif isfield(locs,'chanpos')
        for i = 1 : size(locs.chanpos,1)
            positions{i,1} = num2str(locs.chanpos(i,1),6);
            positions{i,2} = num2str(locs.chanpos(i,2),6);
            positions{i,3} = num2str(locs.chanpos(i,3),6);
            positions{i,4} = ['E_' num2str(i)];
        end
    elseif isfield(locs,'coilpos')
        for i = 1 : size(locs.coilpos,1)
            positions{i,1} = num2str(locs.coilpos(i,1),6);
            positions{i,2} = num2str(locs.coilpos(i,2),6);
            positions{i,3} = num2str(locs.coilpos(i,3),6);
            positions{i,4} = ['E_' num2str(i)];
        end
    end
else
    for i = 1 : size(locs,1)
        positions{i,1} = num2str(locs(i,1),6);
        positions{i,2} = num2str(locs(i,2),6);
        positions{i,3} = num2str(locs(i,3),6);
        positions{i,4} = ['E_' num2str(i)];
    end
end

fid = fopen(filename,'wt');
fprintf(fid,'%s\n',num2str(size(positions,1)));
for i = 1 : size(positions,1)
    fprintf(fid,'%s\t %s\t %s\t %s\t\n',positions{i,:});
end;
fclose(fid);

return