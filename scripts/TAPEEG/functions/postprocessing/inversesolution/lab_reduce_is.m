% Reduce IS-matrix to selected solutionpoints
%
% written by F. Hatz 2012

function ISmatrix = lab_reduce_is(ISmatrix,include)

if ~exist('include','var')
    include = 1:ISmatrix.numsolutionpoints;
end

if isfield(ISmatrix,'x')
    ISmatrix.x = ISmatrix.x(include,:,:);
    ISmatrix.y = ISmatrix.y(include,:,:);
    ISmatrix.z = ISmatrix.z(include,:,:);
end
if isfield(ISmatrix,'matrix')
    ISmatrix.matrix = ISmatrix.matrix(include,:,:);
end
ISmatrix.TSolutionPointName = ISmatrix.TSolutionPointName(include,:);
ISmatrix.numsolutionpoints = length(include);