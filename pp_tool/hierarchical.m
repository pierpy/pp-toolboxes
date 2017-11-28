V=[];


Z = zscore(total_data');

for i = 1:6
    i
    M=[];
    idx = clusterdata(Z', 'linkage', 'ward', 'distance', 'euclidean', 'maxclust', i);
    for k=1:i
        M(:,k)=mean(Z(:,idx==k),2);
        
    end
    V(i)=mean(var(M,[],2));
    
end
count=0;
        kk=0;
        m=1;
        durata=[];
for m=1:4
     for k = 1:size(idx,1)-1
                if idx(k) == m 
                    count=count+1;
                    if idx(k)~=idx(k+1)
                        kk=kk+1;
                        durata(kk)=count;
                        count=0;
                    end
                end
     end
                 meandur(m)=mean(durata);
                 maxdur(m)=max(durata);
                 mindur(m)=min(durata);
end
        