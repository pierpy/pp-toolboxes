function lab_write_montage(Filename,Montage)

if ~exist('Montage','var') | ~isfield(Montage,'chans') | ~isfield(Montage,'label')
    return
end
if isempty(Filename)
    [Montage_file,Montage_filepath] = uiputfile('*.xls;*.xlsx','Select File to store','Montage.xls');
    if Montage_file == 0
        return
    end
    Filename = fullfile(Montage_filepath,Montage_file);
end

for j = 1:size(Montage.chans,1)
    xlsout{j,1} = Montage.label{j,1}; %#ok<AGROW>
    if Montage.chans{j,2} == 1
        xlsout{j,2} = ['AUX' num2str(Montage.chans{j,1})]; %#ok<AGROW>
    else
        xlsout{j,2} = Montage.chans{j,1}; %#ok<AGROW>
    end
    if Montage.chans{j,4} == 1
        xlsout{j,3} = ['AUX' num2str(Montage.chans{j,3})]; %#ok<AGROW>
    else
        xlsout{j,3} = Montage.chans{j,3}; %#ok<AGROW>
    end
end
xlsout = cat(1,{'Montage',Montage.numchans,'channels';Montage.name,'active','reference'},xlsout);
lab_write_xls(Filename,xlsout);
