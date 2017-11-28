function [r, e] = assortativity(g)
% hierarchy     - Assortativity coefficient
%   
%   r = hierarchy(g) find assortative coefficient. r = 1 indicate perfect
%   maxing and r = 0 means there is no assortative maxing. minimun value of
%   r 0 is -1. Random graph tend to negative side.
%
%   [r, e] = hierarchy(g) return assortative maxing matrix
%   

% Ref: M. E. Newman 2003, Maxing pattern in network

r = 1; % perfect maxing
gp = get(g, 'group');
ugp = unique(gp);
ng = length(ugp);

if ng == 1
    return;
end

m = size(g.edges, 1);
e = zeros(ng, ng);
for k = 1:m
    startNodeGroup = find(ugp==g.nodes(g.edges(k,1)).groupid);
    endNodeGroup = find(ugp==g.nodes(g.edges(k,2)).groupid);
    e(startNodeGroup, endNodeGroup) = e(startNodeGroup, endNodeGroup) + 1/m;
end

a = sum(e,1);
b = sum(e,2);

r = trace(e) - a * b / (1 - a * b);