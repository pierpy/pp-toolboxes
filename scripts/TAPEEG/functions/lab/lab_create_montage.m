function montage = lab_create_montage(numchans,header)

montage.name = 'Input structure';
montage.numchans = numchans;
if exist('header','var') & ~isempty(header) & isfield(header,'channels')
    montage.label = cellstr(header.channels);
else
     montage.label = cellstr(num2str((1:numchans)'));
end
for i = 1:numchans
    montage.chans{i,1} = i;
    montage.chans{i,2} = 0;
    montage.chans{i,3} = 0;
    montage.chans{i,4} = 0;
end