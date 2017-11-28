% Write Cartool .spi
%
% filename = lab_write_spi(filename,locs)
%
% written by F. Hatz 2012

function filename = lab_write_spi(filename,locs,fixname)

if ~exist('filename','var')
    [filename,filepath]=uiputfile('*.spi','Select file and path');
    filename = fullfile(filepath,filename);
    clearvars filepath;
end
if ~exist('fixname','var') | fixname ~= 1
    [~,filepath,~,filename] = lab_filename(filename);
    filename = fullfile(filepath,[filename '.spi']);
end

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
            positions{i,4} = ['SP_' num2str(i)];
        end
    elseif isfield(locs,'coilpos')
        for i = 1 : size(locs.coilpos,1)
            positions{i,1} = num2str(locs.coilpos(i,1),6);
            positions{i,2} = num2str(locs.coilpos(i,2),6);
            positions{i,3} = num2str(locs.coilpos(i,3),6);
            positions{i,4} = ['SP_' num2str(i)];
        end
    end
else
    for i = 1 : size(locs,1)
        positions{i,1} = num2str(locs(i,1),6);
        positions{i,2} = num2str(locs(i,2),6);
        positions{i,3} = num2str(locs(i,3),6);
    end
end

fid = fopen(filename,'wt');
if size(positions,2) == 4
    for i = 1 : size(positions,1)
        fprintf(fid,'%s\t%s\t%s\t%s\n',positions{i,:});
    end
else
    for i = 1 : size(positions,1)
        fprintf(fid,'%s\t%s\t%s\n',positions{i,:});
    end
end
fclose(fid);

return