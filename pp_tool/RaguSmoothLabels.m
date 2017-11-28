function [smoothedLabels,bf,org_l] = RaguSmoothLabels(V,M,b,lambda,start,stop)
% function [l,org_l] = RaguSmoothLabels(V,M,b,lambda)
%
% Where  V is the voltage vector(Ne x Nt)
%        M is the microstate vector(Nu x Ne)
%        b is the window size
%       lambda is the non-smoothness penalty factor (should be between 0 and 1)
%

% Copyright 2009-2011 Thomas Koenig
% distributed under the terms of the GNU AFFERO General Public License


if nargin < 5
    start = 1;
end

if nargin < 6
    stop = size(V,2);
end

%if start > 1
%    MSClass(:,:,1:(start-1)) = 0;
%end
%
%if stop < size(MSClass,3)
%    MSClass(:,:,(stop+1):end) = 0;
%end


[Ns,Nt] = size(V);
Nu = size(M,1);

M = M';
fit = V'*M;

if(Nu == 1)
    smoothedLabels = ones(1,Nt);
    bf = fit;
    org_l = smoothedLabels;
    return
end

[mx,l] = max(fit,[],2);
Vvar = sum(V.*V,1);
e = (sum(Vvar) - mx'*mx) / (Nt * (Ns - 1));

if nargout > 3
    org_l = zeros(size(l));
    org_l(start:stop) = l(start:stop);
end

lambda = lambda / (2*b+1);

Oldl = zeros(size(l));
misfit = (repmat(Vvar,Nu,1) - fit'.*fit'.*sign(fit')) / (2*e*(Ns-1));

hits = zeros(Nu,Nt);

while sum(Oldl - l) ~= 0
    Oldl = l;
    for t = (1+b) : (Nt-b)
        for k = 1:Nu
            hits(k,t) = sum(l((t-b):(t+b)) == k);
        end
    end
    penalty = misfit - lambda*hits;
    [mn,l] = min(penalty);
    l = l';
end
smoothedLabels = zeros(size(l));
smoothedLabels(start:stop) = l(start:stop);

if nargout > 1
    bf = zeros(Nt,1);
    for t = start:stop
        bf(t) = fit(t,l(t));
    end
end
