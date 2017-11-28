function MRK = lab_shuffle_markers(MRK,Markers,dosort)
    
    if ~exist('dosort','var')
        dosort = false;
    end
    if ~isfield(MRK,'TYP') | ~isfield(MRK,'POS') | isempty(MRK.POS)
        return
    end
    IDX = [];
    for i = 1:length(Markers)
        tmp = find(strcmp(MRK.TYP,Markers{i}));
        IDX = [IDX tmp(:)']; %#ok<AGROW>
    end
    CLUST.POS = MRK.POS(1,IDX);
    CLUST.DUR = MRK.DUR(1,IDX);
    CLUST.OFF = MRK.OFF(1,IDX);
    CLUST.TYP = MRK.TYP(1,IDX);
    IDX = setdiff(1:length(MRK.POS),IDX);
    MRK.POS = MRK.POS(1,IDX);
    MRK.DUR = MRK.DUR(1,IDX);
    MRK.OFF = MRK.OFF(1,IDX);
    MRK.TYP = MRK.TYP(1,IDX);
    numclust = length(CLUST.POS);
    IDX = randperm(numclust);
    CLUST.TYP = CLUST.TYP(1,IDX);
    CLUST.DUR = CLUST.DUR(1,IDX);
    if dosort == true
        for i = 2:numclust
            CLUST.POS(1,i) = int64(CLUST.POS(1,i-1) + CLUST.DUR(1,i-1));
        end
    end
    MRK = lab_mix_markers(CLUST,MRK);
end