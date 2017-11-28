function [MSClass,MSFit,gfp] = RaguFitMicrostates(gm,MSMaps,Smooth,b,lambda,Start,Stop)

% Copyright 2009-2011 Thomas Koenig
% distributed under the terms of the GNU AFFERO General Public License

MSMaps = NormDimL2(MSMaps,2) / sqrt(size(MSMaps,2));

ng = size(gm,1);
nc = size(gm,2);
ne = size(gm,3);
nt = size(gm,4);

gfp = std(gm,1,3);

MSFit   = zeros(ng,nc,nt);
MSClass = zeros(ng,nc,nt);

for g = 1:ng
    for c = 1:nc
        if Smooth == 1 && size(MSMaps,1) > 1
            [idx,bf] = RaguSmoothLabels(squeeze(gm(g,c,:,:)),MSMaps,b,lambda,Start,Stop );
        else
            [bfall,idxall] = max(MSMaps * squeeze(gm(g,c,:,:)),[],1);
            bf = zeros(1,size(gm,4));
            idx = bf;
            bf(Start:Stop) = bfall(Start:Stop);
            idx(Start:Stop) = idxall(Start:Stop);
        end
        MSClass(g,c,:) = idx;
        MSFit(g,c,:)   = bf;
    end
end
    
MSFit = MSFit / sqrt(ne);
