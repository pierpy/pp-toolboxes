function [ clusteringResults ] = segmentation( self, data, gfp_peaks_indices)
    tic
% CW : Clustering Wrapper
% input : data [ nsamples x nchannels ]
%         minclusters -> min number of clusters to search
%         maxclusters -> max number of clusters to search
%         pmode -> 0: polarity ignored, 1: polarity matters
%         algorithm --> clustering algorithm to perform -> 1: modified k-means (Pascual-Marqui et al., 1995)
%                                                          2: modified AAHC (Murray et al., 2008)
%                                                          3: simple k-means (MatLab function) 
    if nargin < 3
        gfp_peaks_indices = 1;
    end
    if self.confstruct.use_gfp_peaks == 1
        data = data(gfp_peaks_indices,:);
    end
    ClusterNr = self.confstruct.minclusters:self.confstruct.maxclusters;
    numclusters = length(ClusterNr);
    [~, Nchans]=size(data);
    for q = 1 : numclusters
        fprintf(strcat('-cl: ', num2str(q)));
        switch self.confstruct.algorithm
            case 1 % clusteringResultsostates alg.
                 [mps, ind, assigned_cluster_ampl, exp_var] = self.modifiedkmeans(data, ClusterNr(q), 20);
                 clusteringResults.exp_var(q)=exp_var;
            case 2 % AAHC alg.
                 [mps,ind, ~, ~] = self.aahc(data, ClusterNr(q), true);
            otherwise
        end
        if q ~= 1
            % calcolo CV
            [clusteringResults.CV(q)] = self.computecv(data, mps, ind, ClusterNr(q));
            % calcolo W (dispersion measure)
            [clusteringResults.W(q)] = self.computew(data, ind, ClusterNr(q));
            [clusteringResults.m(q)] = self.computem(W(data, ind, ClusterNr(q)), Nchans, ClusterNr(q));
        end
        clusteringResults.template{q} = mps;
        clusteringResults.clusters{q} = ind;
        clusteringResults.gfpmaxima = data;
        clusteringResults.gfp_indices = gfp_peaks_indices;
        if self.confstruct.algorithm == 1
            clusteringResults.amplitude{q} = assigned_cluster_ampl;
        end
    end
    if q ~= 1
        clusteringResults.kl  =  self.computekl(clusteringResults.m);
    end
    toc
end

