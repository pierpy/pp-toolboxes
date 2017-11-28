% Helper script for lab_ICC
%
% Written by F. Hatz 2014

function Result = lab_jackknife(data,func,parameters)

Ntrials = size(data,1);

if exist('parameters','var')
    R0 = func(data,parameters);
else
    R0 = func(data);
end

for i = 1:Ntrials
    include = setdiff(1:Ntrials,i);
    if exist('parameters','var')
        R(i,:) = func(data(include,:),parameters);
    else
        R(i,:) = func(data(include,:));
    end
end

for i = 1:size(R,2)
    tmp = (R(:,i) - R0(1,i)).^2;
    jk(1,i) = (((Ntrials-1) / Ntrials) * sum(tmp));%^0.5;
end

Result.value = R0;
Result.value_low = R0 - 1.96*jk;
Result.value_high = R0 + 1.96*jk;
Result.jackknife = jk;
