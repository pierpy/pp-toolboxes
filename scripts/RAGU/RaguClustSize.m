function y = RaguClustSize(hits)

% Copyright 2009-2011 Thomas Koenig
% distributed under the terms of the GNU AFFERO General Public License

ntime  = size(hits,1);
nruns  = size(hits,2);
y = zeros(1,ntime);

for r = 1:nruns
    cnt = 0;
    for t = 1:ntime
        if hits(t,r)
            cnt = cnt +1;
        else
            if cnt
            y(cnt) = y(cnt) + 1;
                cnt = 0;
            end
        end % end t
        if cnt
            y(cnt) = y(cnt) + 1;
        end
    end
end
