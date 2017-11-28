function [y,y2] = RaguClustSizeV2(hits)
% RaguClustSizeV2
% hits is a timepoint by randomization run matrix
% is the distribution of cluster sizes in the data
% y(1) is the count of clusters with size 1
% y(2) is the count of clusters with size 2
% y(n) is the count of clusters with size n


% Copyright 2009-2011 Thomas Koenig
% distributed under the terms of the GNU AFFERO General Public License

ntime  = size(hits,1); % Number of time points
nruns  = size(hits,2); % Number of randomization run
y = zeros(1,ntime);
y2 = zeros(size(hits));
for r = 1:nruns
    cnt = 0;    % Size of the cluster
    for t = 1:ntime
        if hits(t,r)        % We're in a cluster
            cnt = cnt +1;   % So we increase cluster size by one
        elseif cnt          % We just left a cluster
            y(cnt) = y(cnt) + 1;    % So we increase the counter of the bin with the corresponding cluster number by one
            y2((t-cnt):(t-1),r) = cnt; % Manhattan thing
            cnt = 0;                % and reset the cluster size counter to zero
        end
    end % end t
    if cnt      % Deal with a potential cluster reaching the end of the data
        y(cnt) = y(cnt) + 1;
    end
end