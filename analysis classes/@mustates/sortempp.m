% James Javurek-Humig   (3/15/05)

%INPUTS: 
%   matrix containing eigenvectors in its columns (vectors)
%   matrix containing the corresponding eigenvalues in a diagonal matrix
%   (values)

%The two matricies are sorted from largest to smallest according to the
%eigenvalue matrix.

%RETURNS:
%   matrix containing the sorted eigenvectors
%   diagonal matrix containing the corresponding eigenvalues

%This will give the same result as Matthew Dailey's sortem.m, but is much
%faster for those who are worried about computational time.

function [vectors values] = sortempp(self, vectors, values)

% %this error message is directly from Matthew Dailey's sortem.m
% if nargin ~= 2
%  error('Must specify vector matrix and diag value matrix')
% end;

vals = max(values); %create a row vector containing only the eigenvalues
[svals inds] = sort(vals,'descend'); %sort the row vector and get the indicies
vectors = vectors(:,inds); %sort the vectors according to the indicies from sort
values = max(values(:,inds)); %sort the eigenvalues according to the indicies from sort
values = diag(values); %place the values into a diagonal matrix