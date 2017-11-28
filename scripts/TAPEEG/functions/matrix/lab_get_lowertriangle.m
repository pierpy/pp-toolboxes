function index = lab_get_lowertriangle(nchans,dodiagonal)

if ~exist('dodiagonal','var')
    dodiagonal = true;
end

index = reshape(1:nchans^2,nchans,nchans);
if dodiagonal == true
    index = tril(index,0);
else
    index = tril(index,-1);
end

index = index(index > 0);

end