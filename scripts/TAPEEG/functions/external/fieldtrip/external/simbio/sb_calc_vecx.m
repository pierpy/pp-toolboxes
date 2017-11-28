function vecx = sb_calc_vecx(stiff,vecb,ref);

% SB_CALC_VECX
%
% $Id: sb_calc_vecx.m 8776 2013-11-14 09:04:48Z roboos $

vecdi = zeros(size(stiff,1),1);
vecdi(ref) = 1;
vecva = zeros(size(stiff,1),1);
[stiff, vecb] = sb_set_bndcon(stiff,vecb,vecdi,vecva);
clear vecdi, vecva;
vecx = sb_solve(stiff,vecb);
end
