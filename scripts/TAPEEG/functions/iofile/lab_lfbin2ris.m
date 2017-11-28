function lab_lfbin2ris

[LFbin_file,LFbin_filepath]=uigetfile('*.bin','Select LeadField file');
if isnumeric(LFbin_file)
    return
end
LF = lab_read_LFbin(fullfile(LFbin_filepath,LFbin_file));

RIS.x = permute(LF(:,1,:),[3 1 2]);
RIS.y = permute(LF(:,2,:),[3 1 2]);
RIS.z = permute(LF(:,3,:),[3 1 2]);

[~,~,~,LFbin_fileS] = lab_filename(LFbin_file);
lab_write_ris(fullfile(LFbin_filepath,[LFbin_fileS '.ris']),RIS);