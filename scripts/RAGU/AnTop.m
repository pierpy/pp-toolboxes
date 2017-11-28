function [gfp,p,AllP] = AnTop(data,n_iter,flags)

% Copyright 2009-2011 Thomas Koenig
% distributed under the terms of the GNU AFFERO General Public License

ns = size(data,1);
ne = size(data,2);

quiet = 0;
nrm   = 0;
if nargin > 2
    if ~isempty(findstr(flags,'q'))
        quiet = 1;
    end
    if ~isempty(findstr(flags,'n'))
        nrm = 1;
    end

end

RDiff = zeros(n_iter,1);

if ~quiet
    disp(sprintf('Analyzing %i cases with %i electrodes',ns,ne));
end

if nrm
    data = normr(data); 
end

RDiff(1) = std(mean(data,1),1);
gfp = RDiff(1);
if ~quiet
    h = waitbar(0,'Permuting, please wait...','Position',[100 100 270 60]); % Cosmetics....
end

for i = 2:n_iter                                                             % Here starts the randomization loop
    rdata = zeros(size(data));
    for j = 1:ns
        rdata(j,:) = data(j,randperm(ne));
    end
    
    RDiff(i) = std(mean(rdata,1));
    if ~quiet
        p = sum(RDiff > RDiff(1))/i;                                         % We compute an ad-hoc p for the impatient user
        set(h,'Name',sprintf('ANTOP: %2.0f%%    p: %0.3f',i/n_iter*100,p));    % and display it
        waitbar(i/n_iter,h);    
    end% with some additional cosmetics
end

p = (sum(RDiff >= RDiff(1)))/(n_iter);

if nargout > 2
    pt = 1:n_iter;
    [s,idx] = sort(RDiff,1,'descend');
    AllP(idx) = pt / n_iter;
end

if ~quiet
    close(h);
end