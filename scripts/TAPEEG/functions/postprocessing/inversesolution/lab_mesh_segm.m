% Mesh segments of mri-segmentation (lab_segment_mri)
%
% [bnd,settings] = lab_mesh_segm(Segm,settings,doshort)
%    settings.nvert = number of vertices for mesh
%
% written by F. Hatz 2012

function [bnd,settings] = lab_mesh_segm(Segm,settings,doshort)

disp('    Mesh volumes using fieldtrip code')
if ~exist('doshort','var')
    doshort = false;
end

if ~exist('settings','var') | ~isfield(settings,'nvert')
    settings.nvert = 3000;
end

if isfield(Segm,'transform')
    transform = Segm.transform;
else
    transform = eye(4);
    if isfield(Segm,'originator')
        transform(1:3,4) = -Segm.originator(1,1:3)';
    end
end

if isfield(settings,'mrifile') & ~isempty(settings.mrifile)
    [~,filepath,~,filenameS] = lab_filename(settings.mrifile);
    warning off %#ok<WNOFF>
    mkdir(fullfile(filepath,'Mesh'));
    warning on %#ok<WNON>
    filenameout = fullfile(fullfile(filepath,'Mesh'),filenameS);
    clearvars filepath filenameS
else
    filenameout = [];
end

warning off %#ok<WNOFF>
bnd = [];
nummesh = 1;
tissue = [];
if isfield(Segm,'scalp')
    if length(settings.nvert) >= nummesh
        nvert = settings.nvert(nummesh);
    else
        nvert = settings.nvert(1);
    end
    originator = lab_vol2orig(Segm.scalp);
    [node,face] = triangulate_seg(Segm.scalp,nvert,originator);
    %node = node -repmat(originator,size(node,1),1);
    node = transform *[node ones(size(node,1),1)]';
    node = node(1:3,:)';
    bnd(1,nummesh).pnt = node(:,1:3);
    bnd(1,nummesh).tri = face;
    if ~isempty(filenameout)
        lab_write_spi([filenameout '_scalp.spi'],bnd(1,nummesh).pnt);
    end
    tissue = [tissue cellstr('scalp')];
    nummesh = nummesh+1;
end
if isfield(Segm,'skull') & doshort == false
    if length(settings.nvert) >= nummesh
        nvert = settings.nvert(nummesh);
    else
        nvert = settings.nvert(1);
    end
    originator = lab_vol2orig(Segm.skull);
    [node,face] = triangulate_seg(Segm.skull,nvert,originator);
    %node = node -repmat(originator,size(node,1),1);
    node = transform *[node ones(size(node,1),1)]';
    node = node(1:3,:)';
    bnd(1,nummesh).pnt = node(:,1:3);
    bnd(1,nummesh).tri = face(:,1:3);
    if ~isempty(filenameout)
        lab_write_spi([filenameout '_skull.spi'],bnd(1,nummesh).pnt);
    end
    tissue = [tissue cellstr('skull')];
    nummesh = nummesh+1;
end
if isfield(Segm,'brain') & doshort == false
    if length(settings.nvert) >= nummesh
        nvert = settings.nvert(nummesh);
    else
        nvert = settings.nvert(1);
    end
    originator = lab_vol2orig(Segm.brain);
    [node,face] = triangulate_seg(Segm.brain,nvert,originator);
    %node = node -repmat(originator,size(node,1),1);
    node = transform *[node ones(size(node,1),1)]';
    node = node(1:3,:)';
    bnd(1,nummesh).pnt = node(:,1:3);
    bnd(1,nummesh).tri = face(:,1:3);
    if ~isempty(filenameout)
        lab_write_spi([filenameout '_brain.spi'],bnd(1,nummesh).pnt);
    end
    tissue = [tissue cellstr('brain')];
end
settings.tissue = tissue;
if isfield(settings,'mrifile') & ~isempty(settings.mrifile) & doshort == false
    disp(['      Save result as jpg (' filenameout '_mesh.jpg)'])
    plot.facecolor = [0 1 1;1 0 1;1 1 0;1 0 0;0 1 0;0 0 1];
    plot.filename = [filenameout '_mesh.jpg'];
    plot.plotfaces = true;
    plot.plotedges = false;
    lab_plot_mesh(bnd,plot);
end
warning on %#ok<WNON>

return