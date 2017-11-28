function settings = lab_get_clusteringdata(settings)

if ~exist('settings','var') | isempty(settings)
    settings.criterion = 'distance';
    settings.cutoff = [];
    settings.depth = [];
    settings.distance = 'euclidean';
    settings.linkage = 'average';
    settings.maxclust = [];
    settings.savememory = 'off';
end
    