% Helper script for reading mri-files to set orientation, dimension and
% originator

function mri = lab_set_mricoord(mri,docorrect)

if ~isstruct(mri) | ~isfield(mri,'transform') | ~isfield(mri,'anatomy')
    return
end

if ~exist('docorrect','var')
    docorrect = false;
end

for i = 1:3;
    if max(mri.transform(i,1:3)) >= abs(min(mri.transform(i,1:3))) | min(mri.transform(i,1:3)) > 0
        if find(mri.transform(i,1:3) == max(mri.transform(i,1:3))) == 1
            tmp(i) = 'r';
            tmp2(i) = 1; 
        elseif find(mri.transform(i,1:3) == max(mri.transform(i,1:3))) == 2
            tmp(i) = 'a';
            tmp2(i) = 2; 
        elseif find(mri.transform(i,1:3) == max(mri.transform(i,1:3))) == 3
            tmp(i) = 's';
            tmp2(i) = 3; 
        end
    else
        if find(mri.transform(i,1:3) == min(mri.transform(i,1:3))) == 1
            tmp(i) = 'l';
            tmp2(i) = 1; 
        elseif find(mri.transform(i,1:3) == min(mri.transform(i,1:3))) == 2
            tmp(i) = 'p';
            tmp2(i) = 2; 
        elseif find(mri.transform(i,1:3) == min(mri.transform(i,1:3))) == 3
            tmp(i) = 'i';
            tmp2(i) = 3; 
        end
    end
end
mri.coordold = tmp;
mri.coordsys = 'ras';
[mri.anatomy,T] = affine(mri.anatomy,mri.transform,[1 1 1],0,0,2);
mri.transform = T;
mri.dimensions = [size(mri.anatomy,tmp2(1)) size(mri.anatomy,tmp2(2)) size(mri.anatomy,tmp2(3))];
if docorrect == true
    if abs(mri.transform(1,4)) <= 10 | abs(mri.transform(1,4)) > (mri.dimensions(1)-10)
        mri.transform(1,4) = -round(mri.dimensions(1)/2);
    end
    if abs(mri.transform(2,4)) <= 10 | abs(mri.transform(2,4)) > (mri.dimensions(2)-10)
        mri.transform(2,4) = -round(mri.dimensions(2)/2);
    end
    if abs(mri.transform(3,4)) <= 10 | abs(mri.transform(3,4)) > (mri.dimensions(3)-10)
        mri.transform(3,4) = -round(mri.dimensions(3)/2);
    end
end
mri.originator = -mri.transform(1:3,4)';