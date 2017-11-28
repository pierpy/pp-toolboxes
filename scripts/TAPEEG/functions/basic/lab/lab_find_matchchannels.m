function cfg = lab_find_matchchannels(cfg,calc)

if ~isfield(cfg,'LOCS') | ~isfield(cfg.LOCS,'locs') | isempty(cfg.LOCS.locs) | ...
        ~isfield(cfg.LOCS,'matchlocs') | cfg.LOCS.matchlocs == false
    return
end
    
disp('   Find matching channels of all processed files')

LOCS = cfg.LOCS.locs;
LabelsAll = LOCS.labels;
xlsout = zeros(length(LOCS.labels),size(calc.Filelist,2));
Filenames = cell(1,size(calc.Filelist,2));
for filenr = 1:size(calc.Filelist,2)
    [EEG_file,~,~,Filenames{1,filenr}] = lab_filename(calc.Filelist{1,filenr});
    disp(['   read ' EEG_file ' (only header)'])
    [~,header,cfg] = lab_read_data(calc.Filelist{1,filenr},cfg,true);
    cfg.LOCS.locs = header.locs;
    
    Ilabels = [];
    for i = 1:length(LabelsAll)
        tmp = find(strcmp(header.locs.labels,LabelsAll{i})==1,1);
        if ~isempty(tmp)
            Ilabels = [Ilabels i]; %#ok<AGROW>
        end
    end
    xlsout(Ilabels,filenr) = 1;
end

if isfield(cfg,'settings_path')
    lab_write_els(fullfile(cfg.settings_path,'electrodes-matching.els'),cfg.LOCS.locs);
    xlsout = cat(2,LabelsAll(:),num2cell(xlsout));
    xlsout = cat(1,[{''} Filenames],xlsout);
    lab_write_xls(fullfile(cfg.settings_path,'electrodes-matching.xlsx'),xlsout);
end

end