% data = lab_read_vrb(filename)
% 
% by F. Hatz, Neurology Basel

function data = lab_read_vrb(filename)

delimiter = '\t';
formatSpec = '%s%[^\n\r]';
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter,  'ReturnOnError', false);
fclose(fileID);
data = [dataArray{1:end-1}];