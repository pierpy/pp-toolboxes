function [gd]=generalizedDissimilarityTest(group1, group2, group3)

allConditions=[group1;group2;group3];
meanAllConditions=mean(allConditions);
residual=allConditions-repmat(meanAllConditions,size(allConditions,1),1);
meanRes1=mean(residual(1:size(group1,1),:));
meanRes2=mean(residual(size(group1,1)+1:size(group1,1)+size(group2,1),:));
meanRes3=mean(residual(size(group2,1)+size(group1,1)+1:size(residual,1),:));
meanRes=[meanRes1;meanRes2;meanRes3];

gd=sum(mean(meanRes.^2,2));



