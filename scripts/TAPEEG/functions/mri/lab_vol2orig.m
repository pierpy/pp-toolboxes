% Find originator in mri-file
%
% originator = lab_vol2orig(vol)

function originator = lab_vol2orig(vol)

X = sum(sum(vol,2),3);
Y = sum(sum(vol,1),3);
Z = sum(sum(vol,1),2);

originator(1) = find(X>0,1,'last')/2 + find(X>0,1,'first')/2;
originator(2) = find(Y>0,1,'last')/2 + find(Y>0,1,'first')/2;
originator(3) = find(Z>0,1,'last')/2 + find(Z>0,1,'first')/2;

return