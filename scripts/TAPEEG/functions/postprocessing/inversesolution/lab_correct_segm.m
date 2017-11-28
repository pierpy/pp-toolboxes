% Correct segmented volumes for holes using spm_dilate_erode and fillholes3d
%
% Segm = lab_correct_segm(Segm)
%   Input = output of lab_segment_mri
%
% written by F. Hatz 2012

function Segm = lab_correct_segm(Segm)

disp('    Correct volumes for intersections')

if ~isfield(Segm,'transform')
    return
end
tmp = find(Segm.transform(3,1:3)<0);
if ~isempty(tmp)
    dflip = tmp;
else
    dflip = 0;
    tmp = find(Segm.transform(3,1:3)>0);
end
if tmp == 1
    pmatrix = [2 3 1];
    pmatrix2 = [3 1 2];
elseif tmp == 2
    pmatrix = [1 3 2];
    pmatrix2 = [1 3 2];
else
    pmatrix = [1 2 3];
    pmatrix2 = [1 2 3];
end

if isfield(Segm,'scalp') & ~isfield(Segm,'corrected')
    scalp = Segm.scalp;
    skull = Segm.skull;
    brain = Segm.brain;
    
    % -- create kernels for dilate / erode --
    ekernel3 = cat(3,[0 0 0; 0 1 0; 0 0 0],[0 1 0; 1 1 1; 0 1 0], ...
        [0 0 0; 0 1 0; 0 0 0]);
    % ekernel5 = cat(3,[0 0 0 0 0;0 0 0 0 0;0 0 1 0 0;0 0 0 0 0;0 0 0 0 0], ...
    %     [0 0 0 0 0;0 0 1 0 0;0 1 1 1 0;0 0 1 0 0;0 0 0 0 0], ...
    %     [0 0 1 0 0;0 1 1 1 0;1 1 1 1 1;0 1 1 1 0;0 0 1 0 0], ...
    %     [0 0 0 0 0;0 0 1 0 0;0 1 1 1 0;0 0 1 0 0;0 0 0 0 0], ...
    %     [0 0 0 0 0;0 0 0 0 0;0 0 1 0 0;0 0 0 0 0;0 0 0 0 0]);
    
    % correct for transformation matrix
    
    % -- fill holes --
    scalp = logical(scalp + skull + brain);
    skull = logical(skull + brain);
    [scalp,lowz] = correctvol(scalp,dflip,pmatrix,pmatrix2,-1,3);
    skull = correctvol(skull,dflip,pmatrix,pmatrix2,lowz+5,3);
    brain = correctvol(brain,dflip,pmatrix,pmatrix2,lowz+10,5);
    
    % -- correct by dilate / erode --
    tmp = spm_dilate_erode(double(brain),ekernel3,'dilate');
    skull(tmp==1) = 1;
    tmp = spm_dilate_erode(double(scalp),ekernel3,'erode');
    skull(tmp==0) = 0;
    tmp = spm_dilate_erode(double(skull),ekernel3,'erode');
    brain(tmp==0) = 0;
    
    % -- collect results --
    Segm.scalp = scalp;
    Segm.skull = skull;
    Segm.brain = brain;
    
    % set flag for corrected
    Segm.corrected = 1;
end

if isfield(Segm,'csf') & isfield(Segm,'brain')
    csf = Segm.csf;
    gray = Segm.gray;
    white = Segm.white;
    csf = correctvol(csf,dflip,pmatrix,pmatrix2,10,0);
    gray = correctvol(gray,dflip,pmatrix,pmatrix2,10,0);
    white = correctvol(white,dflip,pmatrix,pmatrix2,10,0);
    csf(Segm.brain==0) = 0;
    gray(Segm.brain==0) = 0;
    white(Segm.brain==0) = 0;
    Segm.csf = csf;
    Segm.gray = gray;
    Segm.white = white;
end

return

function [vol,mask] = correctvol(vol,dflip,pmatrix,pmatrix2,mask,fillholes)

if dflip > 0
    if exist('flipdim') == 2 %#ok<EXIST>
        vol = flipdim(vol,dflip);
    elseif exist('flip') == 2 %#ok<EXIST>
        vol = flip(vol,dflip);
    end
end
vol = permute(vol,pmatrix);

if mask == -1
    mask = find(sum(sum(vol,1),2)> size(vol,1),1)+ 3;
    voltmp = vol(:,:,1:mask);
end
vol(:,:,1:mask) = 1;
if fillholes > 0
    vol = fillholes3d(vol,fillholes);
end
if exist('voltmp','var')
    vol(:,:,1:mask) = voltmp;
    clearvars voltmp
else
    vol(:,:,1:mask) = 0;
end
%vol = smoothbinvol(vol,3);

vol = permute(vol,pmatrix2);
if dflip > 0
    if exist('flipdim') == 2 %#ok<EXIST>
        vol = flipdim(vol,dflip);
    elseif exist('flip') == 2 %#ok<EXIST>
        vol = flip(vol,dflip);
    end
end

return