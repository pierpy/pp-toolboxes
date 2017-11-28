function [b_model, ind, c, gfp] = AAHC(eeg, n_mod, IgnorePolarity)

% Copyright 2009-2011 Thomas Koenig
% distributed under the terms of the GNU AFFERO General Public License


if nargin < 3 
    IgnorePolarity = false;
end

if (size(n_mod,1) ~= 1)
	error('Second argument must be a scalar')
end

%if (size(n_mod,2) ~= 1)
%	error('Second argument must be a scalar')
%end

[n_frame,n_chan] = size(eeg);

% 
h = eye(n_chan)-1/n_chan;
eeg = eeg*h;									% Average reference of data 
gfp=std(eeg,1,2);
Cluster = NormDimL2(eeg,2);

FrameAssignedTo = 1:n_frame;
IsAlive = true(1,n_frame);

GEV = sum(eeg.*Cluster,2);
if (IgnorePolarity == true)
    GEV = abs(GEV);
end

nClusters = sum(IsAlive == true);
tic
while (nClusters) > min(n_mod)
    
    [m,ToRemove] = nanmin(GEV);   % This is the cluster that has to die
    IsAlive(1,ToRemove) = false;    % And we do some bookkeeping
    GEV(ToRemove,1) = NaN;
    
    ToReassign = find(FrameAssignedTo == ToRemove); % These are the orphans
    
    Fit = eeg(ToReassign,:) * Cluster'; % Now this is the fit of the orphans with the remaining clusters
    
    if (IgnorePolarity == true)
        Fit = abs(Fit);
    end
            
    Fit(:,IsAlive == false) = -Inf;
 
    [m,NewAssignment] = max(Fit,[],2);  % We decide where they go
       
    FrameAssignedTo(1,ToReassign) = NewAssignment;  % And do the assignment

    for f = 1:numel(NewAssignment)  % Some cluster center need an update
        ClusterWeDealWith = NewAssignment(f);
        ClusterMembers = FrameAssignedTo == ClusterWeDealWith;
        if (IgnorePolarity == false)
            Cluster(ClusterWeDealWith,:) = NormDimL2(mean(eeg(ClusterMembers,:),1),2);
        else
            
            [pc1,~] = eigs(cov(eeg(ClusterMembers,:)),1);

            
            Cluster(ClusterWeDealWith,:) = NormDimL2(pc1',2);
        end
        NewFit = Cluster(ClusterWeDealWith,:)* eeg(ClusterMembers,:)';

        
        if (IgnorePolarity == true)
            NewFit = abs(NewFit);
        end
        
        GEV(ClusterWeDealWith,1) = sum(NewFit,2);
    end
    
    nClusters = sum(IsAlive == true);

    prc = 1-(nClusters+min(n_mod))/n_frame;

    if numel(n_mod) > 1
        idx = find(n_mod == nClusters);
        if numel(idx) > 0
            b_model{idx} = Cluster(IsAlive,:);
        end
    end
end

if numel(n_mod) == 1
    b_model = Cluster(IsAlive,:);
end
  covm    = eeg * b_model';
  [c,ind] =  max(covm,[],2);
end

function [m,i] = nanmin(v)

idx = find(~isnan(v));
[m,p] = min(v(idx));
i = idx(p);
end

