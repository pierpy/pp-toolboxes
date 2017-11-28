function lab_collect_avgsweeps(searchfolder)
    
disp('Collect AVG sweeps')

if exist('searchfolder','var') & exist(searchfolder,'dir')
    cfg.SEARCH.searchfolder = searchfolder;
end
cfg.SEARCH.searchstring = {'AvgSelect.mrk'};
[cfg,skipprocessing] = lab_set_searchstrings(cfg,1,0,{'strings'});
if skipprocessing == 1;
    return
end
[FileList,cfg] = lab_search_files(cfg);
List = FileList.Filelist;
if isempty(List)
    disp('No files to process found')
    return
else
    List = List(:)';
end

xlsout = {};
progressbar('Collect .mrk-files')
ListN = size(List,2);
for N = 1:ListN
    progressbar(N/ListN);
    Markers = lab_read_mrk(List{N});
    if isfield(Markers,'TYP')
        xlsout{end+1,1} = List{N}(1:end-10); %#ok<AGROW>
        xlsout{end,2} = sum(strcmp(Markers.TYP,'Sweep'));
    end
end

if ~isempty(xlsout)
    lab_write_xls(fullfile(cfg.SEARCH.searchfolder,'AVG_Sweeps_Count.xlsx'),xlsout);
end