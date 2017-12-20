

function [vectors values] = itab_sortempp(vectors, values)

% %this error message is directly from Matthew Dailey's sortem.m
% if nargin ~= 2
%  error('Must specify vector matrix and diag value matrix')
% end;

vals = max(values); %create a row vector containing only the eigenvalues
[svals inds] = sort(vals,'descend'); %sort the row vector and get the indicies
vectors = vectors(:,inds); %sort the vectors according to the indicies from sort
values = max(values(:,inds)); %sort the eigenvalues according to the indicies from sort
values = diag(values); %place the values into a diagonal matrix