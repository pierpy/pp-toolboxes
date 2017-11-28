% Function to correct landmarks
%  To set landmarks in a mri-file, use the function lab_plot_orthoslides.
%  Landmarks are set as 'fixed' or 'correct'. If = 'correct', the position
%  of the landmark is recalculated by this function and the landmark placed
%  between the two indicated other landmarks at the indicated position
%  (default = 50% distance).
%  Landmarks are stored in the file *.lmrk ( = standard Matlab-container-format)
%
% mri = lab_calculate_landmarks(mri)
%
% written by F. Hatz 2014

function mri = lab_calculate_landmarks(mri)

disp('    Correct landmarks')

if ischar(mri)
    [MRI_file,MRI_filepath,~,MRI_fileS] = lab_filename(mri);
    cd(MRI_filepath);
    warning off %#ok<WNOFF>
    mkdir(fullfile(MRI_filepath,'MriSeg'));
    warning off %#ok<WNOFF>
    MRI_file_out = fullfile(fullfile(MRI_filepath,'MriSeg'),MRI_fileS);
    if ~exist([MRI_file_out '_scalp.hdr'],'file')
        cfgtmp = [];
        cfgtmp.mrifile = fullfile(MRI_filepath,MRI_file);
        cfgtmp.SEGcorrect = 1;
        cfgtmp.SEGprobmaps = 0;
        mri = lab_segment_mri(cfgtmp);
        clearvars cfgtmp
    else
        disp(['     read Scalp-Info from ' MRI_file_out '_scalp.hdr'])
        mri = lab_read_mri(fullfile(MRI_filepath,MRI_file));
        mritmp = lab_read_mri([MRI_file_out '_scalp.hdr']);
        mri.scalp = logical(mritmp.anatomy);
        clearvars mritmp
    end
end

if ~isfield(mri,'landmarks') | isempty(mri.landmarks)
    disp('     Abort: no landmarks defined');
    return
end
if ~isfield(mri,'scalp')
    disp('     Abort: missing scalp information');
    return
end

% Find landmarks for wich formula was defined (landmark.calc1
% landmark.calc2 landmark.calcfactor)
Nlm = [];
for i = 1:length(mri.landmarks)
    if strcmp(mri.landmarks(i).mode,'correct')
        Nlm = [Nlm i];
    end
    LandmarkLabel{i,1} = mri.landmarks(i).name;
end

if ~isempty(Nlm)
    % Calculate low resolution scalp mesh
    cfgtmp = [];
    cfgtmp.nvert = 3000;
    [Lbnd,cfgtmp] = lab_mesh_segm(mri,cfgtmp);
    if isfield(cfgtmp,'tissue')
        tmp = find(strcmp(cfgtmp.tissue,'scalp'));
        if isempty(tmp)
            tmp = 1;
        end
    else
        tmp = 1;
    end
    Lbnd = Lbnd(tmp);
    clearvars cfgtmp tmp
end

for N = Nlm
    % intersect scalp-maesh and plane (defined by the three landmarks)
    landmark = mri.landmarks(N);
    pnt1 = mri.landmarks(strcmp(LandmarkLabel,landmark.calc1)).pnt;
    pnt2 = mri.landmarks(strcmp(LandmarkLabel,landmark.calc2)).pnt;
    plane = createPlane(landmark.pnt+1e-10,pnt1+1e-10,pnt2+1e-10);
    disp('     intersect scalp-mesh and plane (please wait)')
    polys = intersectPlaneMesh(plane,Lbnd.pnt,Lbnd.tri);
    polys = polys{1,1};
    clearvars plane
    
    % delete lowest points
    polys = polys(polys(:,3) > (min(polys(:,3))+4),:);
    
    % sort by distance
    polys = cat(1,polys,pnt2,pnt2,pnt1,pnt1);
    distance = lab_distance(polys);
    Idx = size(polys,1);
    tmp = [];
    for i = 1:5
        tmp(end+1,:) = polys(Idx,:);
        distance(Idx,:) = NaN;
        [~,tmp2] = sort(distance(:,Idx));
        tmp2 = tmp2(1:20);
        while length(tmp2) > 1 & (polys(tmp2(1),3) - tmp(end,3)) < 0
            tmp2 = tmp2(2:end);
        end
        if length(tmp2) == 1
            Idx = find(distance(:,Idx) == min(distance(:,Idx)),1);
        else
            Idx = tmp2(1);
        end
        clearvars tmp2
    end
    while size(tmp,1) ~= size(polys,1);
        tmp(end+1,:) = polys(Idx,:);
        distance(Idx,:) = NaN;
        Idx = find(distance(:,Idx) == min(distance(:,Idx)),1);
    end
    polys = tmp;
    clearvars tmp distance Idx i
        
    % calculate line distance
    distance = line_distance(polys);
    tmp = find(distance==0);
    tmp1 = tmp(1:find(diff(tmp)>1,1));
    tmp2 = tmp(find(diff(tmp)>1,1)+1:end);
    if length(tmp1) == 1
        tmp1 = tmp1 + 1;
    else
        tmp1 = tmp1(end);
    end
    if length(tmp2) ~= 1
        tmp2 = tmp2(1);
    else
        tmp2 = tmp2 - 1;
    end
    polys = polys(tmp1+1:tmp2,:);
    clearvars tmp1 tmp2 distance
    
    % Increase resolution to maximal distance < 1mm
    distance = line_distance(polys);
    while max(distance) > 1
        N2 = size(polys,1);
        tmp = zeros(N2-1,3);
        for i = 2:N2
            tmp(i-1,:) = [polys(i,1)+polys(i-1,1) polys(i,2)+polys(i-1,2) polys(i,3)+polys(i-1,3)] ./ 2;
        end
        tmp2 = zeros(2*N2-1,3);
        tmp2(1:2:end,:) = polys;
        tmp2(2:2:end,:) = tmp;
        polys = tmp2;
        clearvars N2 tmp tmp2
        distance = line_distance(polys);
    end
    clearvars tmp tmp2
        
    % calculate pnt
    distance2 = distance;
    for i = 2:length(distance)
        distance2(i) = distance2(i) + distance2(i-1);
    end
    tmp = abs(distance2 - (sum(distance) * landmark.calcfactor / 100));
    tmp = find(tmp == min(tmp),1);
    tmp = polys(tmp+1,:);
    tmp = round(tmp * 100) / 100;
    mri.landmarks(N).pnt = tmp;
    mri.landmarks(N).mode = 'calculated';
    clearvars polys tmp i distance pnt1 pnt2 landmark
end

if exist('MRI_filepath','var')
    % store landmarks
    landmarks = mri.landmarks;
    save(fullfile(MRI_filepath,[MRI_fileS '.lmrk']),'landmarks');
end

end

function distance = line_distance(locs)
   numlocs = size(locs,1);
   distance = zeros(1,numlocs-1);
   for i = 2:numlocs
       distance(1,i-1) = sqrt((locs(i,1)-locs(i-1,1)).^2 + ...
           (locs(i,2)-locs(i-1,2)).^2 + (locs(i,3)-locs(i-1,3)).^2);
   end
end