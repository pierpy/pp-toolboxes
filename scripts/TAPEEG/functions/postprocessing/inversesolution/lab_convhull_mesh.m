% Create convhull-meshes for forward-solution-BEM-calculation
%
% [bndout,settings] = lab_convhull_mesh(bnd,settings)
%
% written by F. Hatz 2013

function [bnd2,settings] = lab_convhull_mesh(bnd,settings)

skipprocessing = 0;

disp('    Hull meshs for BEM')

if ~exist('settings','var') | ~isfield(settings,'SkullDist')
    settings = lab_get_MESH(settings,'BEM');
end

if settings.SkullNr == 0 & isfield(settings,'tissue') & ~isempty(settings.tissue)
    SkullNr = find(strcmp(settings.tissue,'skull'));
else
    SkullNr = settings.SkullNr;
end

if skipprocessing == 0
    nummesh = size(bnd,2);
    thickness = ones(1,nummesh)*settings.OtherDist;
    if SkullNr > 0 & SkullNr <= nummesh
        thickness(SkullNr) = settings.SkullDist;
    end
    if max(bnd(1).pnt(:,1)) > max(bnd(nummesh).pnt(:,1)) | max(bnd(1).pnt(:,3)) > max(bnd(nummesh).pnt(:,3))
        tmp = nummesh:-1:1;
        bnd = bnd(tmp);
        thickness = thickness(tmp);
        if isfield(settings,'tissue')
            settings.tissue = settings.tissue(1,tmp);
        end
        clearvars tmp
    end
    
    disp('      Convhull mesh')
    for i = 1:nummesh
        vertices = bnd(i).pnt;
        numvertices(i) = size(vertices,1);
        kvertices = convhulln(vertices);
        vertices = vertices(kvertices,:);
        vertices = unique(vertices, 'rows');
        faces = projecttri(vertices,'convhull');
        bnd(i).pnt = vertices;
        bnd(i).tri = faces;
    end
    
    disp('      Correct for minimal distances')
    if nummesh > 1
        for i = 2:nummesh
            vertices = bnd(i).pnt;
            faces = bnd(i).tri;
            vertices2 = bnd(i-1).pnt;
            numvertices2 = size(vertices2,1);
            for j = 1:numvertices2
                isgood = 0;
                while isgood == 0
                    tmp = [];
                    ratio = 0.1;
                    while isempty(tmp)
                        [tmp,tmp2,tmp3] = intersectLineMesh3d([0 0 0 vertices2(j,:)*ratio],vertices,faces);
                        tmp = tmp(tmp2>0,:);
                        ratio = 2*ratio;
                        if ratio>1000
                            tmp = 'skip';
                        end
                        %disp(['double ' num2str(j)])
                    end
                    if strcmp(tmp,'skip')
                        rdiff = 0;
                    else
                        rdiff = norm(vertices2(j,:)) + thickness(i) - norm(tmp(1,:));
                    end
                    if rdiff > 0
                        tmp3 = tmp3(tmp2>0);
                        vtmp = vertices(faces(tmp3(1),1),:);
                        vtmpnorm = vtmp ./ norm(vtmp);
                        vertices(faces(tmp3(1),1),:) = vtmp + rdiff*vtmpnorm;
                        vtmp = vertices(faces(tmp3(1),2),:);
                        vtmpnorm = vtmp ./ norm(vtmp);
                        vertices(faces(tmp3(1),2),:) = vtmp + rdiff*vtmpnorm;
                        vtmp = vertices(faces(tmp3(1),3),:);
                        vtmpnorm = vtmp ./ norm(vtmp);
                        vertices(faces(tmp3(1),3),:) = vtmp + rdiff*vtmpnorm;
                    else
                        isgood = 1;
                    end
                end
            end
            bnd(i).pnt = vertices;
        end
    end
    
    disp('      Convhull mesh')
    for i = 1:nummesh
        vertices = bnd(i).pnt;
        kvertices = convhulln(vertices);
        vertices = vertices(kvertices,:);
        vertices = unique(vertices, 'rows');
        faces = projecttri(vertices,'convhull');
        bnd(i).pnt = vertices;
        bnd(i).tri = faces;
    end
    
    disp('      Re-tesselate mesh using brainstrom code')
    for i = 1:nummesh
        [bnd(i).pnt,bnd(i).tri] = tess_refine(bnd(i).pnt,bnd(i).tri,3);
        [bnd(i).pnt,bnd(i).tri] = tess_refine(bnd(i).pnt,bnd(i).tri,3);
        [bnd(i).pnt,bnd(i).tri] = tess_remesh(double(bnd(i).pnt),numvertices(i));
        %     bnd(i).pnt = round(bnd(i).pnt*10)/10;
    end
    
    bnd2 = bnd;
    if isfield(settings,'fileout') & ~isempty(settings.fileout)
        [~,filepath,~,filenameS] = lab_filename(settings.fileout);
        warning off %#ok<WNOFF>
        mkdir(fullfile(filepath,'Mesh'));
        warning on %#ok<WNON>
        filenameout = fullfile(fullfile(filepath,'Mesh'),filenameS);
        clearvars filepath filenameS
        disp(['      Save result as jpg (' filenameout '_BEMmesh.jpg)'])
        plot.facecolor = [0 1 1;1 0 1;1 1 0;1 0 0;0 1 0;0 0 1];
        plot.plotedges = false;
        plot.plotfaces = true;
        plot.filename = [filenameout '_BEMmesh.jpg'];
        lab_plot_mesh(bnd2,plot);
        clearvars plot
        for i = 1:size(bnd2,2)
            if isfield(settings,'tissue')
                lab_write_spi([filenameout '_BEM_' settings.tissue{i} '.spi'],bnd2(i).pnt);
            else
                lab_write_spi([filenameout '_BEM_' num2str(i) '.spi'],bnd2(i).pnt);
            end
        end
    end
else
    bnd2 = bnd;
end