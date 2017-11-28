% Read Cartool .spi
%
% LOCS = lab_read_spi(filename)
%
% written by F. Hatz 2012

function LOCS = lab_read_spi(filename)

if ~exist('filename','var') | ~exist(filename,'file')
    [SPI_file,SPI_filepath]=uigetfile({'*.spi'},'Select spi file');
    filename = fullfile(SPI_filepath,SPI_file);
end

delimiter = '\t';
formatSpec = '%f%f%f%s%[^\n\r]';
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'EmptyValue' ,NaN, 'ReturnOnError', false);
fclose(fileID);
clearvars filename delimiter formatSpec fileID ans;

LOCS.x = dataArray{:,1}';
LOCS.y = dataArray{:,2}';
LOCS.z = dataArray{:,3}';
LOCS.labels = dataArray{:,4}';
LOCS = lab_locs2sph(LOCS);

return