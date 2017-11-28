% Write Cartool .rois
%
% lab_write_rois(filename,ROIS)
%
% written by F. Hatz 2012

function lab_write_rois(filename,ROIS)
[~,filepath,~,filename] = lab_filename(filename);
filename = fullfile(filepath,[filename '.rois']);

fid = fopen(filename,'wt');
fprintf(fid,[ROIS.version '\n']);
fprintf(fid,[num2str(ROIS.numsolutionpts) '\n']);
fprintf(fid,[num2str(ROIS.numrois) '\n']);
for i = 1 : ROIS.numrois
    fprintf(fid,[ROIS.labels{i} '\n']);
    fprintf(fid,[num2str(ROIS.solutionpts{i}) '\n']);
end;
fclose(fid);