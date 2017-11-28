function b_model = RaguEEG_Mod_AAHC(eeg,n_mod,ProgBar, IgnorePolarity, max_n)

% Copyright 2009-2011 Thomas Koenig
% distributed under the terms of the GNU AFFERO General Public License

if nargin < 3
    ProgBar = true;
end

if nargin < 4 
    IgnorePolarity = false;
end

if (size(n_mod,1) ~= 1)
	error('Second argument must be a scalar')
end

%if (size(n_mod,2) ~= 1)
%	error('Second argument must be a scalar')
%end

[n_frame,n_chan] = size(eeg);

if nargin > 4
    idx = randperm(n_frame);
    eeg = eeg(idx(1:max_n),:);
    [n_frame,n_chan] = size(eeg);
end

h = eye(n_chan)-1/n_chan;
eeg = eeg*h;									% Average reference of data 

if (ProgBar == true)
    hndl = waitbar(0,sprintf('Fitting %i-%i microstates (AAHC), please wait...',n_mod(1),n_mod(end)));
end

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
    
    if (rem(nClusters,10) == 0) && (ProgBar == true)
        waitbar(prc,hndl);
        set(hndl,'Name',sprintf('Remaining time: %01.0f:%02.0f min',floor(toc()*(1/prc)/60),rem(toc()*(1/prc-1),60)));
    end
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

if (ProgBar == true)
    close(hndl);
end
end

function [m,i] = nanmin(v)

idx = find(~isnan(v));
[m,p] = min(v(idx));
i = idx(p);
end

