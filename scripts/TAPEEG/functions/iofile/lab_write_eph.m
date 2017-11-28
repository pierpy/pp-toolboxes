% Write Cartool .eph
%
% lab_write_eph(filename,data,header)
%
% written by F. Hatz 2012

function lab_write_eph(filename,data,header)
[~,filepath,~,filenameS] = lab_filename(filename);
filename = [filenameS '.eph'];

head=[size(data,1) size(data,2) header.samplingrate];
dlmwrite(fullfile(filepath,filename),head,'\t');
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

% Write Marker-file
if isfield(header,'events')
    Marker_file=fullfile(filepath,[filename '.mrk']);
    lab_write_mrk(Marker_file,header);
end

% Write individual freqbands
if isfield(header,'IFREQ') & ~isempty(header.IFREQ)
    IFREQ = header.IFREQ; %#ok<NASGU>
    save(fullfile(filepath,[filenameS '.ifreq']),'IFREQ');
end

return