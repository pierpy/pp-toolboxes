function ll = Ragu_ComputeLoreta(map,lor)

if size(map,1) ~= lor.nchan
    error('Channel number mismatch in function Ragu_ComputeLoreta');
end

l = zeros(3,size(map,2),size(lor.spinv,3));

l(1,:,:) = map'* squeeze(lor.spinv(:,1,:));
l(2,:,:) = map'* squeeze(lor.spinv(:,2,:));
l(3,:,:) = map'* squeeze(lor.spinv(:,3,:));

ll = sum(l.*l,1);
ll = shiftdim(ll,1)';

