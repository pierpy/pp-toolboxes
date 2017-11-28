% Reduce ROIS-information to selected solutionpoints
%
% rois = lab_reduce_rois(rois,selection)
%
% written by F. Hatz 2012

function rois = lab_reduce_rois(rois,selection)

if size(selection,2) == 1
    selection = selection';
end

if max(selection) <= rois.numsolutionpts
    rois.solutionptsAll = [];
    for i = 1:size(rois.solutionpts,2)
        tmp = zeros(1,rois.numsolutionpts);
        tmp(rois.solutionpts{1,i}) = 1;
        tmp = tmp(1,selection);
        rois.solutionpts{1,i} = find(tmp==1);
        rois.solutionptsAll = union(rois.solutionptsAll,rois.solutionpts{1,i});
    end
    if size(rois.solutionptsAll,1) > 1
        rois.solutionptsAll = rois.solutionptsAll';
    end
    rois.numsolutionpts = size(selection,2);
end