% Write digitizer information of Elekta system to .els
%
% lab_write_digitizer(header,filename)
%
% written by F. Hatz 2012

function lab_write_digitizer(filename,header)
[~,filepath,~,filename] = lab_filename(filename);
filename = fullfile(filepath,[filename '.spi']);

for i = 1:size(header.orig.info.dig,2)
    positions{i,1} = header.orig.info.dig(1,i).r(1,1)*1000;
    positions{i,2} = header.orig.info.dig(1,i).r(2,1)*1000;
    positions{i,3} = header.orig.info.dig(1,i).r(3,1)*1000;
    positions{i,4} = [num2str(header.orig.info.dig(1,i).kind) '-' num2str(header.orig.info.dig(1,i).ident)];
end

fid = fopen(filename,'wt');
for i = 1 : size(positions,1)
    fprintf(fid,'%10.4f\t %10.4f\t %10.4f\t %s\t\n',positions{i,:});
end;
fclose(fid);