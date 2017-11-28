function lab_write_mapping(Filename,mappings)

if ~exist('mappings','var') | ~isfield(mappings,'mappings') | ~isfield(mappings,'mappingstitle')
    return
end
if isempty(Filename)
    [Mappings_file,Mappings_filepath] = uiputfile('*.xls;*.xlsx','Select File to store','Mappings.xls');
    if Mappings_file == 0
        return
    end
    Filename = fullfile(Mappings_filepath,Mappings_file);
end

xlsout = [mappings.mappingstitle mappings.mappings'];

if isfield(mappings,'mappingsChannelsFile')
    allchans = ones(1,mappings.mappingsChannelsFile);
elseif isfield(mappings,'mappingsChannels')
    allchans = ones(1,mappings.mappingsChannels);
else
    maxchan = 0;
    for i = 1:size(mappings.mappings,2)
        maxchan = max(maxchan,max(mappings.mappings{1,i}));
    end
    allchans = ones(1,maxchan);
end
for i = 1:size(mappings.mappings,2)
    allchans(1,mappings.mappings{1,i}) = 0;
    xlsout{i,2} = num2str(xlsout{i,2});
end
xlsout{end+1,1} = 'no region';
xlsout{end,2} = num2str(find(allchans == 1));
lab_write_xls(Filename,xlsout);
