% Extract lower triangle of matrix
%
% matrix = lab_tril2matrix(connections,Nchans)
%
% Written by F. Hatz 2013

function matrix = lab_tril2matrix(connections,Nchans)

if ~exist('Nchans','var')
    Nchans = floor((size(connections,1)*2)^0.5);
elseif ~isnumeric(Nchans) | length(Nchans) > 1
    tmp = length(Nchans);
    Nchans = tmp;
    clearvars tmp
end
if isempty(connections)
    matrix = [];
    return
end
if size(connections,1) == 1
    connections = connections';
end
if mod(size(connections,1),Nchans^2/2 + Nchans/2) == 0
    dodiag = true;
    Nconn = Nchans^2/2 + Nchans/2;
elseif mod(size(connections,1),Nchans^2/2 - Nchans/2) == 0
    dodiag = false;
    Nconn = Nchans^2/2 - Nchans/2;
else
    matrix = [];
    return
end
connections = reshape(connections,[Nconn size(connections,1)/Nconn*size(connections,2)]);
Nmatrix = size(connections,2);
if dodiag == true
    idx = tril(true(Nchans,Nchans));
else
    idx = tril(true(Nchans,Nchans),-1);
end
matrix = zeros(Nchans,Nchans,Nmatrix);
for i = 1:Nmatrix
    Mtmp = zeros(Nchans,Nchans);
    Mtmp(idx) = connections((i-1)*Nconn+1:i*Nconn);
    tmp = diag(Mtmp);
    Mtmp = Mtmp + Mtmp';
    Mtmp(1:Nchans+1:end) = tmp;
    matrix(:,:,i) = Mtmp;
end