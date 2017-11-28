function [w] = W(data,clust_ind,ClusterNr)
     % input arguments
     % data = the input data (number of time instances * number of channels)
     % clust_ind = vector with assignement of each time frame with the most correlated template (number of time instances * 1)
     % ClusterNr = number of current clusters

     % output arguments
     % cv for the current cluster
     Clusterstmp = clust_ind;
     for j = 1:ClusterNr
        D(j)= sum(sum(triu(squareform(pdist(data(Clusterstmp==j,:),'euclidean')))))/(2*size(data(Clusterstmp==j,:),2));
     end
     w=sum(D); 
end