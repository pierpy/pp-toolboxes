function [cv] = itab_computecv(templates, clust_ind, ClusterNr)
     % input arguments
     % data = the input data (number of time instances * number of channels)
     % templates = templates calculated for the current cluster (number of templates * number of channels)
     % clust_ind = vector with assignement of each time frame with the most correlated template (number of time instances * 1)
     % ClusterNr = number of current clusters

     % output arguments
     % cv for the current cluster
    [Ntr, Nchans]=size(data);
    cv_coeff = (Nchans-1)/(Nchans-1-ClusterNr);
    cv_coeff = cv_coeff^2;
    templateTemp = templates;
    Clusterstmp = clust_ind;
    tmp = zeros(size(data,1),1);
    for j = 1:size(data,1)
        tmp(j,1) = data(j,:)*data(j,:)' - (templateTemp(Clusterstmp(j,1),:)*data(j,:)').^2;
    end
    cv = (sum(tmp)/(Ntr*(Nchans-1)))*cv_coeff;
end