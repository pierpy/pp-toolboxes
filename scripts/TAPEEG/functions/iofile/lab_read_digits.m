function digits = lab_read_digits(filename)

if ~exist('filename','var')
    [filename,filepath] = uigetfile('*.fif;*.fiff','select MEG data file');
    filename = fullfile(filepath,filename);
    clearvars filepath
end

[hdr] = fiff_setup_read_raw(filename);
tmp = hdr.info.dig;
digits = zeros(size(tmp,2),3);
for i = 1:size(tmp,2)
   digits(i,:) = tmp(1,i).r';
end
digits = digits * 1000;