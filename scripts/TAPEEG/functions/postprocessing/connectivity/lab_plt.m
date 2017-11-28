% Calculate PLT (phase locking time)
%
% plt = lab_plt(phases,samplingrate)
%
% written by F.Hatz 2014

function plt = lab_plt(phases,samplingrate)

nstep = ceil((size(phases,2)-1)/4095);

SumDiff = zeros(size(phases,1),size(phases,1));
for i = 1:nstep
    Sstart = (i-1)*4095+1;
    if i~=nstep
        Sstop = i*4095+1;
        if Sstop > size(phases,2)
            Sstop = size(phases,2);
        end
    else
        Sstop = size(phases,2);
    end
    phasetmp = phases(:,Sstart:Sstop);
    diffs = repmat(permute(phasetmp,[1 3 2]),[1 size(phases,1) 1]) - repmat(permute(phasetmp,[3 1 2]),[size(phases,1) 1 1]);
    diffs = sign(sin(diffs));
    diffs(diffs == 0) = 1;
    diffs = sum(diff(diffs,[],3)~=0,3);
    SumDiff = SumDiff+diffs;
end
plt = 1 - exp(-(SumDiff*samplingrate).^-1 * size(phases,2));
plt(1:size(plt,1)+1:end) = 0;
