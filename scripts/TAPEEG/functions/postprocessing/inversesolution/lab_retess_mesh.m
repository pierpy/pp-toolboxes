% Re-tesselate Mesh using a sphere and project the shpere to the volume
%
% bnd = lab_retess_mesh(bnd,numvertices)

function bnd = lab_retess_mesh(bnd,numvertices)

disp('    Re-tess mesh')

nummesh = size(bnd,2);
for i = 1:nummesh
    vertices = bnd(i).pnt;
    faces = bnd(i).tri;
    [pnt2, tri2] = mysphere(numvertices); % this is a regular triangulation
    for j = 1:size(pnt2,1)
        [tmp tmp2] = intersectLineMesh3d([0 0 0 pnt2(j,:)],vertices,faces);
        tmp = round(tmp);
        tmp = tmp(tmp2>0,:);
        ratio = 0.1;
        while isempty(tmp)
            [tmp,tmp2] = intersectLineMesh3d([0 0 0 pnt2(j,:)*ratio],vertices,faces);
            tmp = tmp(tmp2>0,:);
            ratio = 2*ratio;
            if ratio>1000
                tmp = 'skip';
            end
            %disp(['double ' num2str(j)])
        end
        if ~strcmp(tmp,'skip')
            pnt2(j,:) = tmp(1,:);
        end
    end
    bnd(i).pnt = pnt2;
    bnd(i).tri = tri2;
end