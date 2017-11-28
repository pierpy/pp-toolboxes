% Calculate PLV (phase locking value)
%
% plv = lab_plv(data,method)
%
%    data = matrix (chans x numtf)
%
% written by F. Hatz 2013

function plv = lab_plv(data,method)

if ~exist('method','var')
    method = 3;
end

hilb = hilbert(data')';

if method == 1
    %method1
    hilbdiff = zeros(size(data,1),size(data,1),size(data,2));
    for j = 1:size(hilb,2)
        hilbdiff(:,:,j) = hilb(:,j) * ((hilb(:,j).').^-1);
    end
    plv = abs(sum((hilbdiff ./ abs(hilbdiff)),3)) ./ size(hilbdiff,3);
elseif method == 2
    % method2 (fastest!!)
    hilbdiff = repmat(permute(angle(hilb),[1 3 2]),[1 size(hilb,1) 1]) - repmat(permute(angle(hilb),[3 1 2]),[size(hilb,1) 1 1]);
    plv = abs(sum(exp(1i*hilbdiff),3))/size(hilbdiff,3);
else
    %method3
    plv = zeros(size(data,1),size(data,1));
    for i = 1:size(data,1)
        for j = 1:size(data,1)
            if i ~= j
                hilbdiff = angle(hilb(i,:)) - angle(hilb(j,:));
                plv(i,j) = abs(sum(exp(1i*hilbdiff),2))/size(hilbdiff,2);
            else
                plv(i,j) = 0;
            end
        end
    end
end