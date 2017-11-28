function [p]=groupAnova(group1, group2, group3, group4, group5, group6)

allConditions=[group1;group2;group3;group4;group5;group6];
meanAllConditions=mean(allConditions);
residual=allConditions-repmat(meanAllConditions,size(allConditions,1),1);
meanRes1=mean(residual(1:size(group1,1),:));
meanRes2=mean(residual(size(group1,1)+1:size(group1,1)+size(group2,1),:));
meanRes3=mean(residual(size(group2,1)+size(group1,1)+1:size(residual,1),:));
meanRes=[meanRes1;meanRes2;meanRes3];
gdReal=sum(mean(meanRes.^2,2));

for k=1:5000
    ind=randperm(size(allConditions,1));
    allConditions_perm=allConditions(ind,:);
    meanAllConditions_perm=mean(allConditions_perm);
    residual_perm=allConditions_perm-repmat(meanAllConditions_perm,size(allConditions_perm,1),1);
    
    meanRes1_perm=mean(residual_perm(1:size(group1,1),:));
    meanRes2_perm=mean(residual_perm(size(group1,1)+1:size(group1,1)+size(group2,1),:));
    meanRes3_perm=mean(residual_perm(size(group2,1)+size(group1,1)+1:size(residual_perm,1),:));
    
    meanRes_perm=[meanRes1_perm;meanRes2_perm;meanRes3_perm];
    gdPerm(k)=sum(mean(meanRes_perm.^2,2));
    clear ind meanCperm meanLperm
end

a=find(gdPerm>=gdReal);

p=size(a,2)/5000;




