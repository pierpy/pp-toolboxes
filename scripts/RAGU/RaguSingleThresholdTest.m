function [threshold,p]= RaguSingleThresholdTest(rd) 

% Copyright 2009-2011 Thomas Koenig
% distributed under the terms of the GNU AFFERO General Public License

n_it = size(rd.PTanova,4);

%m = squeeze(min(rd.PTanova,[],3));
m = squeeze(mean(rd.PTanova,3));
 
ms = sort(m,3,'ascend');
 
idx = floor(n_it * rd.Threshold);
if (isempty(idx)||(idx < 1))
    threshold = 0;
   p = 1;
   return;
end
    
threshold = squeeze(ms(:,:,idx));

p = sum(m <= repmat(m(:,:,1),[1,1,n_it]),3) / n_it;
 
 
 