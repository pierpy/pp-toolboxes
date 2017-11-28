function [template_to_return]=mean_templates(data, n_clusters, zscore)
 %   template_to_return: mean subjects wise templates calculated with the
 %                       algoritm described in Koenig et al., 1999
 %   data:              (n_subjects*n_clusters)*n_channels matrix.
 %   n_clusters:        number of clusters (microstates).
 %   zscore:            0 if no z-score; 1 instead.
 
 if zscore
    Z = zscore(data);
 else
    Z=data;
 end
    
for rep = 1:20
    
    idx = randperm(size(Z,1), n_clusters);
    template = Z(idx,:);
    covm1    = Z * template';
    [c,ind] =  max(covm1,[],2);
     for k=1:n_clusters
        ii = find(ind==k);
        template(k,:)=mean(Z(ii,:));
     end
     covm2 = Z*template';
     [new_c,new_ind] =  max(covm2,[],2);
     length(find(ind ~= new_ind))

    for iter = 1:50
            disp(strcat('iterazione: ', num2str(iter)))
            for k=1:n_clusters
                ii2 = find(new_ind==k);
                template_new(k,:)=mean(Z(ii2,:));
            end
            covm = Z*template_new';
            [new_c,new_ind] =  max(covm,[],2);
            stability(iter)=length(find(ind ~= new_ind));
            
    end
    TEMPLATES{rep}=template_new;
    V(rep)=mean(var(template_new,[],2));
end  
    minVar = find(V==min(V));
    template_to_return = TEMPLATES{minVar};
end

