function [ NMicr ] = CW( data, minclusters, maxclusters, pmode, algorithm )
% CW : Clustering Wrapper
% input : data [ nsamples x nchannels ]
%         minclusters -> min number of clusters to search
%         maxclusters -> max number of clusters to search
%         pmode -> 0: polarity ignored, 1: polarity matters
%         algorithm --> clustering algorithm to perform -> 1: modified k-means (Pascual-Marqui et al., 1995)
%                                                          2: modified AAHC (Murray et al., 2008)
%                                                          3: simple k-means (MatLab function) 

ClusterNr = minclusters:maxclusters;
numclusters = length(ClusterNr);
[Ntr, Nchans]=size(data);


for q=1:numclusters
    %disp(strcat('-cl: ', num2str(q)));

    switch algorithm
        case 1 % Nmicrostates alg.
            %disp('modified version of k-mea  ns')
             [mps, ind, assigned_cluster_ampl, exp_var] = NMicrostates(data, ClusterNr(q), pmode, 20);
             NMicr.exp_var(q)=exp_var;
        case 2 % AAHC alg.
            %disp('modified version of AAHC')
             [mps,ind, c, gfp] = AAHC(data,ClusterNr(q),true);
%                           NMicr.exp_var(q) = (sum(c*gfp'))/(sum(gfp.^2));
        case '3' % simple k-means alg
            %disp('simple k-means')
             [ind, mps]=kmeans(data,q);
        otherwise
            %disp('other value')
    end

    if numclusters ~= 1
        % calcolo CV
        [NMicr.CV(q)] = CV(data,mps,ind,ClusterNr(q));
        % calcolo W (dispersion measure)
        [NMicr.W(q)] = W(data,ind,ClusterNr(q));
        [NMicr.m(q)] = M(W(data,ind,ClusterNr(q)),Nchans,ClusterNr(q));
    end

    NMicr.template{q}=mps;
    NMicr.clusters{q}=ind;
    NMicr.originaldata = data;
    NMicr.gfpmaxima = data;
%     NMicr.gfp_maxima_allData = allData;
        % NMicr.amplitude{q}=assigned_cluster_ampl;
end
if numclusters ~= 1
    NMicr.kl  =  KL(NMicr.m);
end

end

