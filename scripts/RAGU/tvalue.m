function t = tvalue(data1,data2,dim,vartype)

if nargin < 3
    dim = 1;
end

if nargin < 4
    vartype = 'e';
end


if isempty(data2)
    t = mean(data1,dim) ./ std(data1,0,dim) * sqrt(size(data1,dim));
else
    switch vartype
        case 'e'
            v = (var(data1,0,dim) * (size(data1,dim)-1) +  var(data2,0,dim) * (size(data2,dim)-1)) / (size(data1,dim) +  size(data2,dim) -2);  % Correct
%            v = (var(data1,0,dim) * (size(data1,dim)-1) +  var(data2,0,dim) * (size(data2,dim)-1)) / (size(data1,dim) +  size(data2,dim));  % Roberto
            t = (mean(data1,dim) -mean(data2,dim))./ sqrt(v * (1/size(data1,dim) + 1/ size(data2,dim)));
        case 'u'
            t = (mean(data1,dim) -mean(data2,dim))./ sqrt(var(data1,0,dim)./size(data1,dim) + var(data2,0,dim)./size(data2,dim));
        otherwise
            error('tvalue: Variance type falsly defined');
    end
end