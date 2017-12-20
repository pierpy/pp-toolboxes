function [ clusteringResults ] = itab_segmentation(data, gfp_peaks_indices)
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
    if  confstruct.use_gfp_peaks == 1
        data = data(gfp_peaks_indices,:);
    end
    ClusterNr =  confstruct.minclusters: confstruct.maxclusters;
    numclusters = length(ClusterNr);
    [~, Nchans]=size(data);
    for q = 1 : numclusters
        fprintf(strcat('-cl: ', num2str(q)));
        switch  confstruct.algorithm
            case 1 % clusteringResultsostates alg.
                 [mps, ind, assigned_cluster_ampl, exp_var] =  itab_modifiedkmeans(data, ClusterNr(q), 20);
                 clusteringResults.exp_var(q)=exp_var;
            case 2 % AAHC alg.
                 [mps,ind, ~, ~] =  itab_aahc(data, ClusterNr(q), true);
            otherwise
        end
        if q ~= 1
            % calcolo CV
            [clusteringResults.CV(q)] =  itab_computecv(data, mps, ind, ClusterNr(q));
            % calcolo W (dispersion measure)
            [clusteringResults.W(q)] =  itab_computew(data, ind, ClusterNr(q));
            [clusteringResults.m(q)] =  itab_computem( computew(data, ind, ClusterNr(q)), Nchans, ClusterNr(q));
        end
        clusteringResults.template{q} = mps;
        clusteringResults.clusters{q} = ind;
        clusteringResults.gfpmaxima = data;
        clusteringResults.gfp_indices = gfp_peaks_indices;
        if  confstruct.algorithm == 1
            clusteringResults.amplitude{q} = assigned_cluster_ampl;
        end
    end
    if q ~= 1
        clusteringResults.kl  =   itab_computekl(clusteringResults.m);
    end
    toc
end

