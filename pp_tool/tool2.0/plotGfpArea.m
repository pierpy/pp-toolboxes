



x = 1:size(gfp,2);        
y = gfp;  
baseLine = 0.2;  
index = 30:40
index1 = find(cluster == 1);
index2 = find(cluster == 2);
index3 = find(cluster == 3);
index4 = find(cluster == 4);

for m = 1 : 4
    
    ind = find(cluster == m);
    for i = 1:size(ind) 
        index = ind(i);
        plot(x,y,'b'); 
        hold on;                                     
        h1 = fill(x(index([1 1:end end])),[baseLine y(index) baseLine],'r','EdgeColor','none');
    end
    
end

