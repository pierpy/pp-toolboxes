function lab_shuffle_clusters
    
set.SEARCH.searchstring = {'*.mrk'};
Files = lab_search_files(set);
Filelist = Files.Filelist;
if isempty(Filelist)
    return
end

for filenr = 1:length(Filelist)
    [~,Filepath,Format,FilenameS] = lab_filename(Filelist{filenr});
    if ~strcmp(Format,'mrk');
        break
    end
    disp(['Read and Shuffle ' FilenameS]);
    MRK = lab_read_mrk(Filelist{filenr});
    Markers = unique(MRK.TYP);
    IDX = [];
    for i = 1:length(Markers)
        if strfind(Markers{i},'Clust')
            IDX = [IDX i]; %#ok<AGROW>
        end
    end
    Markers = Markers(IDX);
    MRK = lab_shuffle_markers(MRK,Markers);
    header.events = MRK;
    lab_write_mrk(fullfile(Filepath,[FilenameS '_shuffle.mrk']),header);
end