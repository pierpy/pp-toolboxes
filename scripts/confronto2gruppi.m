clear all
close all
clc

meanC=mean(D_controlli);
meanL=mean(D_right);
diff=meanC-meanL;
GFP_diff=computegfp(diff);

MM=[D_controlli;D_right];


for k=1:5000
ind=randperm(size(MM,1));
MM_perm=MM(ind,:);
meanCperm=mean(MM_perm(1:20,:));
meanLperm=mean(MM_perm(21:end,:));
diffPerm=meanCperm-meanLperm;
GFP_diffPerm(k)=computegfp(diffPerm);
clear ind meanCperm meanLperm
end

a=find(GFP_diffPerm>=GFP_diff);

p=size(a,2)/5000