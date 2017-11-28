function [DATA,PLOT] = lab_plot_sort_data(DATA,PLOT)
    
Flag = false(size(PLOT,1),1);    
for i = 1:size(PLOT,1)
    if size(PLOT,2) > 1 & i > 1 & strcmp(PLOT(i-1,1).Mode,'Connections') & ~strcmp(PLOT(i,1).Mode,'Connections')
        Flag(i,1) = true;
        for j = 1:size(PLOT,2)
            PLOT(i,j).AddPlot = true;
        end
    end
end

Flag2 = find(Flag == 0);
Flag2 = Flag2-1;
Flag2 = setdiff(Flag2,0);
Flag2 = union(Flag2,length(Flag));
Idx = [];
Nr = 1;
for i = 1:length(Flag2)
    tmp = Nr:(Flag2(i)*size(PLOT,2));
    tmp = reshape(tmp,length(tmp)/size(PLOT,2),size(PLOT,2));
    Idx = cat(1,Idx,tmp);
    Nr = max(tmp(:))+1;
end
Idx = Idx(:);
PLOT(Idx) = PLOT;
PLOT = PLOT(:);
DATA(Idx) = DATA;
DATA = DATA(:);