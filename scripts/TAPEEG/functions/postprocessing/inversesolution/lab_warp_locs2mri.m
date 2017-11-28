% Warp electrodes to mri with defined landmarks (landmark-names = electrode-labels)
%    for every electrodes a transformation is generated using the three
%    closest landmark-electrodes (resize, translation, rotation, shear)
%
% [mri,LOCS,settings] = lab_warp_locs2mri(mri,LOCS,scalp,settings)
%
% written by F. Hatz 2014

function [mri,LOCS,settings] = lab_warp_locs2mri(mri,LOCS,scalp,settings)

if ~exist('settings','var')
    settings = [];
end
if ~exist('scalp','var')
    scalp = [];
end
if ~exist('LOCS','var')
    LOCS = [];
end
if ~exist('mri','var')
    mri = [];
end
if ~isfield(mri,'anatomy') | ~isfield(LOCS,'x')
    return
end

[mri,settings] = lab_mri_landmarks(mri,LOCS,settings);
if ~isfield(mri,'landmarks') | isempty(mri.landmarks)
    LOCS = [];
    return
else
    pause(0.2);
end

if isempty(scalp)
    % Calculate high resolution mesh
    cfgtmp = [];
    cfgtmp.nvert = 20000;
    [bnd,cfgtmp] = lab_mesh_segm(mri,cfgtmp);
    if isfield(cfgtmp,'tissue')
        tmp = find(strcmp(cfgtmp.tissue,'scalp'));
        if isempty(tmp)
            tmp = 1;
        end
    else
        tmp = 1;
    end
    scalp = bnd(tmp).pnt;
    scalptri = bnd(tmp).tri;
    clearvars cfgtmp bnd tmp
else
    scalptri = lab_pnt2faces(scalp);
end

% Calculate landmarks for wich formula was defined (landmark.calc1
% landmark.calc2 landmark.calcfactor)
mri = lab_calculate_landmarks(mri);

% find landmarks in LOCS
LandmarkIndex = zeros(size(LOCS.x,2),5);
for i = 1:length(mri.landmarks)
    tmp = find(strcmp(LOCS.labels,mri.landmarks(i).name));
    if ~isempty(tmp)
        LandmarkIndex(tmp(1),1) = 1;
        LandmarkIndex(tmp(1),2) = i;
        pnt = mri.landmarks(i).pnt;
        distance = lab_distance(scalp,pnt);
        LandmarkIndex(tmp(1),3:5) = scalp(find(distance==min(distance),1),:);
    end
end
clearvars pnt distance
if sum(LandmarkIndex(:,1)) < 3
    LOCS = [];
    disp('   Abort: not enough landmarks defined or labels of LOC-file not matching')
    return
end

% calculate distance 2d
disp('    Warp Locs to MRI')
LOCS = lab_locs2sph(LOCS);
[x,y] = pol2cart(LOCS.theta,LOCS.radius);
locs = [x' y' ones(size(x,2),1)];
distance = lab_distance(locs);
distance(:,LandmarkIndex(:,1)==0) = Inf;
locs = [LOCS.x' LOCS.y' LOCS.z'];
locs2 = locs;
Result.store = [];
Result.locs = {};
for i = 1:size(locs,1);
    if LandmarkIndex(i,1) == 1
        locs2(i,:) = LandmarkIndex(i,3:5);
    else
        [~,tmp] = sort(distance(i,:));
        tmp = tmp(1:3);
        tmp = sort(tmp);
        if ~isempty(Result.store)
            tmp2 = [];
            for j = 1:size(Result.store,2)
                if isequal(tmp',Result.store(:,j))
                    tmp2 = j;
                end
            end
        else
            tmp2 = [];
        end
        if ~isempty(tmp2)
            locstmp = Result.locs{tmp2(1)};
        else
            locstmp = calculate_newlocs(locs,tmp,LandmarkIndex(tmp,3:5));
            Result.store(:,end+1) = tmp';
            Result.locs{end+1} = locstmp;
        end
        distance2 = lab_distance(scalp,locstmp(i,:));
        tmp = find(distance2==min(distance2),1);
        if distance2(tmp) < 3
            locs2(i,:) = scalp(find(distance2==min(distance2),1),:);
        else
            PntProj = intersectLineMesh3d([0 0 0 locstmp(i,:)*2],scalp,scalptri);
            distance3 = lab_distance(PntProj,locstmp(i,:));
            locs2(i,:) = PntProj(find(distance3==min(distance3),1),:);
            clearvars PntProj distance3
        end
        clearvars tmp tmp2 locstmp distance2
    end
end
LOCS.x = locs2(:,1)';
LOCS.y = locs2(:,2)';
LOCS.z = locs2(:,3)';

end

function distance = line_distance(locs)
   numlocs = size(locs,1);
   distance = zeros(1,numlocs-1);
   for i = 2:numlocs
       distance(1,i-1) = sqrt((locs(i,1)-locs(i-1,1)).^2 + ...
           (locs(i,2)-locs(i-1,2)).^2 + (locs(i,3)-locs(i-1,3)).^2);
   end
end

function locsout = calculate_newlocs(locs,index,PntT)
   PntO = locs(index,:);
   angleO = anglePoints3d(PntO(1,:),PntO(2,:),PntO(3,:));
   angleT = anglePoints3d(PntT(1,:),PntT(2,:),PntT(3,:));
   distO1 = line_distance(PntO(1:2,:));
   distO2 = line_distance(PntO(2:3,:));
   distT1 = line_distance(PntT(1:2,:));
   distT2 = line_distance(PntT(2:3,:));
   
   PntO2 = [distO1 0 0 ; 0 0 0 ; cos(angleO)*distO2 sin(angleO)*distO2 0];
   [~,~,T] = procrustes(PntO2,PntO,'Scaling',false,'Reflection',false);
   Transform = eye(4);
   Transform(1:3,1:3) = T.T';
   Transform(1:3,4) = mean(T.c,1);
   Transform2 = eye(4);
   Transform2(1,1) = distT1 / distO1;
   Transform2(2,2) = (sin(angleT)*distT2) / (sin(angleO)*distO2);
   Transform2(3,3) = ((distT1 / distO1) + (distT2 / distO2)) / 2;
   YDist = sin(angleT) * distT2;
   XShift = (cos(angleO)*distO2)*(distT1 / distO1) - cos(angleT)*distT2;
   
   locs2 = (Transform2 * Transform * ([locs ones(size(locs,1),1)])')';
   for i = 1:size(locs2,1)
       if locs2(i,2) >= YDist
           locs2(i,1) = locs2(i,1) - XShift;
       else
           locs2(i,1) = locs2(i,1) - (XShift * (locs2(i,2) / YDist));
       end
   end
   PntT2 = locs2(index,1:3);
   % angleT2 = anglePoints3d(PntT2(1,:),PntT2(2,:),PntT2(3,:));
   % distT21 = line_distance(PntT2(1:2,:));
   % distT22 = line_distance(PntT2(2:3,:));
   
   [~,~,T] = procrustes(PntT,PntT2,'Scaling',false,'Reflection',false);
   Transform = eye(4);
   Transform(1:3,1:3) = T.T';
   Transform(1:3,4) = mean(T.c,1);
   locs2 = Transform * locs2';
   locsout = locs2(1:3,:)';
end
