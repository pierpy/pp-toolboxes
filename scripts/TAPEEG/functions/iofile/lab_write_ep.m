% Write Cartool .ep
%
% lab_write_ep(filename,data,header)
%
% written by F. Hatz 2012

function lab_write_ep(filename,data,header)
[~,filepath,~,filenameS] = lab_filename(filename);
filename = [filenameS '.ep'];

dlmwrite(fullfile(filepath,filename),data','delimiter','\t','-append');

if exist('header','var') & isfield(header,'ref_chan')
    % Write EEGinfo-file (*.txt)
    lab_write_eeginfo(fullfile(filepath,filename),header)
end
if exist('header','var') & isfield(header,'locs')    
    % Write loc file
    ELS_file = fullfile(filepath,[filenameS '.els']);
    ELS_file = lab_write_locs(ELS_file,header,'bad');
    if ~isempty(ELS_file)
        % Write *.LM-file
        fidout=fopen(fullfile(filepath,[filenameS '.lm']),'w');
        fprintf(fidout,[filename native2unicode([13 10])]);
        fprintf(fidout,[ELS_file native2unicode([13 10])]);
        fclose(fidout);
    end
end

return