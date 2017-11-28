% Read Cartool .spi
%
% locs = lab_read_spi(filename)
%
% written by F. Hatz 2012

function locs = lab_read_pts(filename)

if ~exist('filename','var') | ~exist(filename,'file')
    [SPI_file,SPI_filepath]=uigetfile({'*.spi;*.hdr;*.nii'},'Select spi file');
    filename = fullfile(SPI_filepath,SPI_file);
end

delimiter = '\t';
formatSpec = '%f%f%f%s%[^\n\r]';
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'EmptyValue' ,NaN, 'ReturnOnError', false);
fclose(fileID);
clearvars filename delimiter formatSpec fileID ans;

locs.x = -dataArray{:,1}';
locs.y = dataArray{:,2}';
locs.z = dataArray{:,3}';
locs.labels = dataArray{:,4}';
locs = lab_locs2sph(locs);

return