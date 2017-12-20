
function [ self ] = createmsstr (confstruct)
% funzione che legge tutta la cartella dei file e crea una cell
% con tutti i soggetti e tutte le condizoni
self.confstruct = confstruct;
if ~ isfield(self.confstruct, 'substract_column_at_start')
    self.confstruct.substract_column_at_start = 1;
end
if ~ isfield(self.confstruct, 'method_GFPeak')
    self.confstruct.method_GFPeak = 'GFPL2';
end
if ~ isfield(self.confstruct, 'use_gfp_peaks')
    self.confstruct.use_gfp_peaks = 1;
end
if ~ isfield(self.confstruct, 'setall1')
    self.confstruct.setall1 = 1;
end
if ~ isfield(self.confstruct, 'normalize')
    self.confstruct.normalize = 1;
end
if ~ isfield(self.confstruct, 'minclusters')
    self.confstruct.minclusters = 1;
end
if ~ isfield(self.confstruct, 'maxclusters')
    self.confstruct.maxclusters = 15;
end
if ~ isfield(self.confstruct, 'pmode')
    self.confstruct.pmode = 0;
end
if ~ isfield(self.confstruct, 'algorithm')
    self.confstruct.algorithm = 1;
end
if ~ isfield(self.confstruct, 'similarity_measure')
    self.confstruct.similarity_measure = 'correlation';
end
if ~ isfield(self.confstruct, 'maxima_method')
    self.confstruct.maxima_method = 'simple';
end
if ~ isfield(self.confstruct, 'debuG')
    self.confstruct.debug = 0;
end
if ~ isfield(self.confstruct, 'minDistance')
    self.confstruct.minDistance=0;
end
if ~ isfield(self.confstruct, 'minAmplitude')
    self.confstruct.minAmplitude=0;
end
if ~ isfield(self.confstruct, 'showpeaks')
    self.confstruct.showpeaks=0;
end
if ~ isfield(self.confstruct, 'downsample')
    self.confstruct.downsample = 0;
end
if  ~ isfield(self.confstruct, 'segment_min_len')
    self.confstruct.segment_min_len = 3;
end
if ~ isfield(self.confstruct, 'tf')
    error('tf is a required field')
end
if  ~ isfield(self.confstruct, 'fs')
    error('fs is a required field')
end

end




