% Write electrodes localications in .sfp format
%
% lab_write_sfp(locs,filename)
%
% written by F. Hatz 2012

function lab_write_sfp(filename,locs)
[~,filepath,~,filename] = lab_filename(filename);
filename = fullfile(filepath,[filename '.sfp']);

if isstruct(locs)
    for i = 1:size(locs.x,2)
        positions{i,1} = locs.labels{1,i};
        positions{i,2} = num2str(-locs.y(1,i),6);
        positions{i,3} = num2str(locs.x(1,i),6);
        positions{i,4} = num2str(locs.z(1,i),6);
    end
else
    for i = 1 : size(locs,1)
        positions{i,1} = ['E_' num2str(i)];
        positions{i,2} = num2str(-locs(i,2),6);
        positions{i,3} = num2str(locs(i,1),6);
        positions{i,4} = num2str(locs(i,3),6);
    end
end

fid = fopen(filename,'wt');
for i = 1 : size(positions,1)
    fprintf(fid,'%s %s %s %s\n',positions{i,:});
end;
fclose(fid);

return