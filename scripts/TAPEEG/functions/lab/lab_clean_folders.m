function lab_clean_folders(searchfolder)

if ~exist('searchfolder','var')
    searchfolder = uigetdir('Select folder');
end

Filelist = lab_search(searchfolder,'._',true);
for i = 1:length(Filelist)
    Filename = lab_filename(Filelist{i});
    if strcmp(Filename(1:2),'._')
        delete(Filelist{i})
    end
end