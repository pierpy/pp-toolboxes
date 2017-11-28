% 3D elimination of borderzones in atlas file
%
% atlas = lab_remove_borderzones(atlas)
%
% written by F. Hatz 2012

function atlas = lab_remove_borderzones(atlas)

fprintf('    Remove borders of regions: ')
Rvol = zeros(size(atlas.anatomy));
nreg = unique(atlas.anatomy(:));
nreg = nreg(nreg>0);
for i = 1:length(nreg)
    vol = zeros(size(atlas.anatomy));
    vol(atlas.anatomy == nreg(i)) = 1;
    vol(atlas.anatomy == 0) = 1;
    vol2 = zeros(size(atlas.anatomy));
    spm_smooth(vol,vol2,[1 1 1]);
    vol2(vol2<1) = 0;
    vol2(atlas.anatomy == 0) = 0;
    if sum(vol2(:)==1) < 100
        fprintf('-')
        vol2 = zeros(size(atlas.anatomy));
        spm_smooth(vol,vol2,[0.5 0.5 0.5]);
        vol2(vol2<1) = 0;
        vol2(atlas.anatomy == 0) = 0;
    end
    if sum(vol2(:)==1) < 100
        fprintf('-')
        vol2 = vol;
        vol2(atlas.anatomy == 0) = 0;
    end
    Rvol(vol2==1) = 1;
    clearvars vol2 vol
    fprintf('.')
end
atlas.anatomy(Rvol==0) = 0;
disp(':')

end