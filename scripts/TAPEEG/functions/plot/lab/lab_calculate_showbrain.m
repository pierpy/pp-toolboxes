% Helper file for lab_plot_IS
%
% script generates from a mri-file a template brain for source plotting, 
% additionally a mri-atlas-file is needed
%
% Written by F. Hatz 2013

function Brain = lab_calculate_showbrain(Brain,dodiag)

store = 0;
if ~exist('dodiag','var')
    dodiag = false;
end
if ~exist('Brain','var')
    Brain = [];
end
if ischar(Brain)
    Mode = Brain;
    Brain = [];
else
    Mode = '';
end
if ~isfield(Brain,'Template')
    if isempty(Mode)
        Mode = questdlg('Template-file or MRI-file','Input','MRI-file','Template-file','Template-file');
    end
    if strcmp(Mode,'Template-file')
        [Brain_file,Brain_filepath]=uigetfile('*.mat','Select MNIbrain.mat');
        if ~isempty(Brain_file) & Brain_file ~= 0
            cd(Brain_filepath);
            Brain = load(fullfile(Brain_filepath,Brain_file));
            Brain.Brain_file = fullfile(Brain_filepath,Brain_file);
        else
            Brain = [];
            return
        end
    else
        Brain.Template = [];
        Prompt = {'MRI-file','mri';'Atlas-file','atlas'};
        Formats = [];
        Formats.type = 'edit';
        Formats.format = 'file';
        Formats.items = {'*.hdr;*.nii'};
        Formats.size = 300;
        Formats(end+1,1).type = 'none';
        Formats(end+1,1).type = 'edit';
        Formats(end,1).format = 'file';
        Formats(end,1).items = {'*.hdr;*.nii'};
        Formats(end,1).size = 300;
        [Brain.Template,Cancelled] = inputsdlg(Prompt,'Select MRI / Atlas',Formats,Brain.Template);
        if Cancelled == 1 | ~exist(Brain.Template.mri,'file') | ~exist(Brain.Template.atlas,'file')
            disp('   no valid files selected')
            Brain = [];
            return
        end
        Brain.Brain_file = Brain.Template.mri;
    end
end
if isfield(Brain.Template,'atlas') & ischar(Brain.Template.atlas)
    progressbar('Match Atlas to MRI');
    Brain.Template.atlas = lab_match_atlas2mri(Brain.Template.atlas,Brain.Template.mri);
end
if ~isfield(Brain.Template,'labels')
    if isfield(Brain.Template.atlas,'labels') & ~isempty(Brain.Template.atlas.labels)
        Brain.Template.labels = Brain.Template.atlas.labels;
    elseif size(unique(Brain.Template.atlas.anatomy),1) == 117
        Brain.Template.labels = lab_get_aal;
    else
        [labels_file,labels_filepath]=uigetfile('*.xls;*.xlsx','Select xls-file with labels');
        if ~isempty(labels_file) & labels_file ~= 0
            cd(labels_filepath);
            if ispc
                [~,Brain.Template.labels] = xlsread(fullfile(labels_filepath,labels_file));
            else
                [~,Brain.Template.labels] = xlsread(fullfile(labels_filepath,labels_file),1,'','basic');
            end
            Brain.Template.labels(:);
        else
            tmp = setdiff(unique(Brain.Template.atlas.anatomy),0);
            for i = 1:length(tmp)
                Brain.Template.labels{i,1} = ['Region_' num2str(tmp(i))];
            end
        end
    end
end

if ~isfield(Brain,'faces') & isfield(Brain.Template,'mri') & isfield(Brain.Template,'atlas') & isfield(Brain.Template,'labels')
    if size(unique(Brain.Template.atlas.anatomy),1) ~= (size(Brain.Template.labels,1) + 1)
        tmp = setdiff(unique(Brain.Template.atlas.anatomy),0);
        Brain.Template.labels = cellstr(num2str(tmp(:)));
        clearvars tmp
    end
    
    % Segment MRI
    progressbar('Segment MRI');
    if ischar(Brain.Template.mri)
        disp('   Segment MRI')
        settings.mrifile = Brain.Template.mri;
        settings.SEGbrain = 1;
        settings.SEGprobmaps = 0;
        settings.SEGcorrect = 1;
        Brain.Template.mri = lab_segment_mri(settings);
        clearvars settings
    elseif ~isfield(Brain.Template.mri,'brain')
        Brain.Template.mri = lab_segment_mri(Brain.Template.mri);
    end
    
    % Set originator
    if isfield(Brain.Template.mri,'originator')
        Brain.Template.originator = Brain.Template.mri.originator;
    elseif isfield(Brain.Template.mri,'transform')
        Brain.Template.originator = abs(Brain.Template.mri.transform(1:3,4)');
    else
        Brain.Template.originator = [floor(size(Brain.Template.mri.brain,1)/2) floor(size(Brain.Template.mri.brain,2)/2) floor(size(Brain.Template.mri.brain,3)/2)];
    end
    
    % Mesh brain
    progressbar('Create brain mesh');
    disp('   Create brain mesh')
    tmp = zeros(size(Brain.Template.mri.gray));
    tmp(Brain.Template.mri.gray) = 1;
    tmp(Brain.Template.mri.white) = 1;
    brain = zeros(size(tmp,1)+1,size(tmp,2),size(tmp,3));
    brain(1:Brain.Template.originator(1),:,:) = tmp(1:Brain.Template.originator(1),:,:);
    brain(end-Brain.Template.originator(1)+1:end,:,:) = tmp(end-Brain.Template.originator(1)+1:end,:,:);
    ekernel3 = cat(3,[0 0 0; 0 1 0; 0 0 0],[0 1 0; 1 1 1; 0 1 0],[0 0 0; 0 1 0; 0 0 0]);
    tmp = spm_dilate_erode(double(brain),ekernel3,'erode');
    tmp = spm_dilate_erode(double(tmp),ekernel3,'dilate');
    brain(tmp==0) = 0;
    CC = bwconncomp(brain);
    CC = CC.PixelIdxList';
    if size(CC,1) > 2
        for i = 3:size(CC,1)
            brain(CC{i}) = 0;
        end
    end
    [faces,vertices]=isosurface(brain);
    [faces,vertices] = reducepatch(faces,vertices,0.1);
    vertices = [vertices(:,2) vertices(:,1) vertices(:,3)];
    %bnd.tri = faces;
    %bnd.pnt = vertices;
    %lab_plot_mesh(bnd);
    Brain.vertices = vertices - repmat(Brain.Template.originator,[size(vertices,1) 1]);
    Brain.faces = faces;
    Brain.originator = [0 0 0];
    clearvars  CC i tmp ekernel3 vertices faces
    
    % Match vertices and regions
    progressbar('Match vertices and atlas-regions');
    disp('   Match vertices and regions')
    tmp =  inputdlg('Error range (voxels)','Error range',[1 25],{'4'});
    if ~isempty(tmp) & ~isempty(tmp{1,1});
        Brain.Template.errorrange = round(str2num(tmp{1,1})); %#ok<ST2NM>
    else
        disp('    no valid error range selected, set to 4')
        Brain.Template.errorrange = 4;
    end
    tmp = zeros(size(Brain.Template.atlas.anatomy));
    tmp2 = setdiff(unique(Brain.Template.atlas.anatomy),0);
    for i = 1:length(tmp2)
        tmp(Brain.Template.atlas.anatomy == tmp2(i)) = i;
    end
    clearvars tmp2
    atlas = zeros(size(tmp,1)+1,size(tmp,2),size(tmp,3));
    atlas(1:Brain.Template.originator(1),:,:) = tmp(1:Brain.Template.originator(1),:,:);
    atlas(end-Brain.Template.originator(1)+1:end,:,:) = tmp(end-Brain.Template.originator(1)+1:end,:,:);
    Numdots = ceil(size(Brain.vertices,1)/100);
    for i = 1:size(Brain.vertices,1)
        tmp = floor(Brain.vertices(i,:)) + Brain.Template.originator;
        tmp2 = atlas(tmp(1),tmp(2),tmp(3));
        step = 1;
        while tmp2 == 0 & step <= Brain.Template.errorrange
            tmp = floor(Brain.vertices(i,:)) + Brain.Template.originator;
            if tmp(1)+step > size(atlas,1)
                tmp(1) = size(atlas,1)-step;
            end
            if tmp(2)+step > size(atlas,2)
                tmp(2) = size(atlas,2)-step;
            end
            if tmp(3)+step > size(atlas,3)
                tmp(3) = size(atlas,3)-step;
            end
            if tmp(1)-step < 1
                tmp(1) = 1+step;
            end
            if tmp(2)-step < 1
                tmp(2) = 1+step;
            end
            if tmp(3)-step < 1
                tmp(3) = 1+step;
            end
            tmp2 = atlas(tmp(1)-step:tmp(1)+step,tmp(2)-step:tmp(2)+step,tmp(3)-step:tmp(3)+step);
            tmp3 = setdiff(unique(tmp2(:)),0);
            maxtmp = [];
            for j = 1:length(tmp3)
                if isempty(maxtmp) | length(find(tmp2(:)==tmp3(j))) > maxtmp
                    maxtmp = length(find(tmp2(:)==tmp3(j)));
                    maxtmp2 = tmp3(j);
                end
            end
            if ~isempty(maxtmp)
                tmp2 = maxtmp2;
            else
                tmp2 = 0;
            end
            step = step + 1;
        end
        Brain.Template.vertex(i,1) = tmp2;
        if mod(i,Numdots) == 0
            progressbar(i/size(Brain.vertices,1));
        end
    end
    clearvars tmp tmp2 tmp3 step i j
    
    % Find center coordinates of regions
    progressbar('Find center coordinates of regions');
    disp('   Find center coordinates of regions')
    tmp = setdiff(unique(atlas),0);
    for i = 1:size(Brain.Template.labels,1)
        [x,y,z] = ind2sub(size(atlas), find(atlas==tmp(i)));
        x2 = mean(x);
        y2 = mean(y);
        z2 = mean(z);
        for j = 1:size(x,1)
            distance(1,j) = ((x(j,1) - x2)^2 + (y(j,1) - y2)^2 + (z(j,1) - z2)^2)^0.5; %#ok<AGROW>
        end
        tmp2 = find(distance == min(distance),1);
        Brain.Template.xyz(i,:) = [x(tmp2),y(tmp2),z(tmp2)] - Brain.Template.originator;
        clearvars tmp2 j x y z distance x2 y2 z2
    end
    clearvars i tmp
    
    % Create meshs for regions
    progressbar('Find center coordinates of regions');
    fprintf('   Create meshs for regions')
    ekernel3 = cat(3,[0 0 0; 0 1 0; 0 0 0],[0 1 0; 1 1 1; 0 1 0],[0 0 0; 0 1 0; 0 0 0]);
    brain2 = spm_dilate_erode(double(brain),ekernel3,'erode');
    tmp = setdiff(unique(atlas),0);
    for i = 1:length(tmp)
        Vregion = zeros(size(atlas));
        Vregion(atlas==tmp(i)) = 1;
        Vregion = spm_dilate_erode(double(Vregion),ekernel3,'dilate');
        Vregion = spm_dilate_erode(double(Vregion),ekernel3,'erode');
        Vregion(brain2==0) = 0;
        CC = bwconncomp(Vregion);
        CC = CC.PixelIdxList';
        if size(CC,1) > 1
            if length(CC{1}) > length(CC{2})*1.5
                for j = 2:size(CC,1)
                    Vregion(CC{j}) = 0;
                end
                clearvars j
            elseif size(CC,1) > 2
                for j = 3:size(CC,1)
                    Vregion(CC{j}) = 0;
                end
                clearvars j
            end
        end
        [Rfaces,Rvertices]=isosurface(Vregion,0.5);
        if ~isempty(Rvertices)
            [Rfaces,Rvertices] = reducepatch(Rfaces,Rvertices,0.1);
            Rvertices = [Rvertices(:,2) Rvertices(:,1) Rvertices(:,3)];
        end
        Brain.Template.regions(1,i).faces = Rfaces;
        Brain.Template.regions(1,i).vertices = Rvertices - repmat(Brain.Template.originator,[size(Rvertices,1) 1]);
        clearvars Rfaces Rvertices CC Vregion
        fprintf('.')
        progressbar(i/length(tmp));
    end
    disp(':')
    clearvars brain2 ekernel3 Vregion tmp i
    progressbar(1);
end

if ~isfield(Brain,'verticesL')
    progressbar('Separate left and right side of brain');
    disp('   Separate left and right side of brain')
    if isfield(Brain,'originator')
        tmp = 1;
        orig = Brain.originator(1)+1;
    else
        tmp(1,1) = sum(Brain.vertices(:,1));
        tmp(1,2) = sum(Brain.vertices(:,2));
        tmp(1,3) = sum(Brain.vertices(:,3));
        tmp = abs(tmp);
        if min(tmp) > 100
            orig = (max(Brain.vertices(:,1)) - min(Brain.vertices(:,1)))/2 + min(Brain.vertices(:,1));
            tmp = 1;
        else
            tmp = find(tmp == min(tmp));
            orig = mean(Brain.vertices(:,tmp));
        end
    end
    Brain.nodeL = find(Brain.vertices(:,tmp)<= orig);
    Brain.nodeR = find(Brain.vertices(:,tmp)>= orig);
    Brain.verticesL = Brain.vertices(Brain.nodeL,:);
    Brain.verticesR = Brain.vertices(Brain.nodeR,:);
    Brain.facesL = [];
    Brain.facesR = [];
    if isfield(Brain,'Template') & isfield(Brain.Template,'regions')
        for i = 1:size(Brain.Template.regions,2)
            if ~isempty(Brain.Template.regions(1,i).vertices)
                if max(Brain.Template.regions(1,i).vertices(:,tmp)) <= orig+3
                    Brain.Template.regionsLR{i} = 'L';
                elseif min(Brain.Template.regions(1,i).vertices(:,tmp)) >= orig-3
                    Brain.Template.regionsLR{i} = 'R';
                else
                    Brain.Template.regionsLR{i} = 'R';
                end
            else
                Brain.Template.regionsLR{i} = 'none';
            end
        end
    end
    Numdots = ceil(size(Brain.faces,1)/100);
    for i = 1:size(Brain.faces,1)
        if intersect(Brain.faces(i,1),Brain.nodeL) & intersect(Brain.faces(i,2),Brain.nodeL) & intersect(Brain.faces(i,3),Brain.nodeL)
            Brain.locL(i,1) = true;
            Brain.locR(i,1) = false;
            Brain.facesL(i,1) = find(Brain.nodeL == Brain.faces(i,1));
            Brain.facesL(i,2) = find(Brain.nodeL == Brain.faces(i,2));
            Brain.facesL(i,3) = find(Brain.nodeL == Brain.faces(i,3));
        else
            Brain.locL(i,1) = false;
        end
        if intersect(Brain.faces(i,1),Brain.nodeR) & intersect(Brain.faces(i,2),Brain.nodeR) & intersect(Brain.faces(i,3),Brain.nodeR)
            Brain.locR(i,1) = true;
            Brain.locL(i,1) = false;
            Brain.facesR(i,1) = find(Brain.nodeR == Brain.faces(i,1));
            Brain.facesR(i,2) = find(Brain.nodeR == Brain.faces(i,2));
            Brain.facesR(i,3) = find(Brain.nodeR == Brain.faces(i,3));
        else
            Brain.locR(i,1) = false;
        end
        if mod(i,Numdots) == 0
            progressbar(i/size(Brain.faces,1));
        end
    end
    if ~isempty(Brain.facesL)
        Brain.facesL = Brain.facesL(Brain.locL,:);
    end
    if ~isempty(Brain.facesR)
        Brain.facesR = Brain.facesR(Brain.locR,:);
    end
    clearvars tmp i
    
    store = 1;
end

if ~isfield(Brain,'labels') | dodiag == true
    if ~isfield(Brain,'labels') | isempty(Brain.labels)
        Brain.labels = Brain.Template.labels;
        Labels = {};
    else
        Labels = Brain.labels;
    end
    LabelsList = setdiff(Brain.Template.labels,Labels);
    LabelsList = cat(1,Labels(:),LabelsList(:));
    
    Prompt = cell(0,2);
    Formats = [];
    
    Prompt(end+1,:) = {'Import regions',''};
    Formats(end+1,1).type = 'button';
    Formats(end,1).style = 'pushbutton';
    Formats(end,1).size = [100 25];
    Formats(end,1).callback = {@load_labels,'labels','labels','@labels'};
    
    Prompt(end+1,:) = {'Save regions',''};
    Formats(end+1,1).type = 'button';
    Formats(end,1).style = 'pushbutton';
    Formats(end,1).size = [100 25];
    Formats(end,1).callback = {@save_labels,[],'labels'};
    
    Prompt(end+1,:) = {'Regions','labels'};
    Formats(end+1,1).type = 'list';
    Formats(end,1).style = 'listbox';
    Formats(end,1).format = 'input';
    Formats(end,1).items = LabelsList;
    Formats(end,1).limits = [0 4];
    Formats(end,1).size = [160 450];
    Formats(end,1).span = [10 2];
    
    [Brain,Cancelled] = inputsdlg(Prompt,'Select regions',Formats,Brain,2);
    if Cancelled == 1
        Brain = [];
        return
    elseif length(Brain.labels) == 1
        disp('Error: Please select more than one region')
        Brain = [];
        return
    else
        docalc = false;
        [tmp1,tmp2,tmp3] = intersect(Labels,Brain.labels);
        if length(tmp1) ~= length(Brain.labels)
            docalc = true;
        elseif max(abs(tmp2-tmp3)) ~= 0
            docalc = true;
        end
        if docalc == true
            Brain.labels = Brain.labels(:);
            progressbar('Define selected atlas-regions');
            Brain.Template2labels = zeros(size(Brain.labels,1),1);
            for i = 1:size(Brain.labels,1)
                tmp = find(~cellfun('isempty',strfind(Brain.Template.labels,strtrim(Brain.labels{i,1}))));
                if ~isempty(tmp)
                    Brain.Template2labels(i,1) = tmp;
                end
            end
            for i = 1:size(Brain.Template.vertex,1)
                tmp = find(Brain.Template2labels == Brain.Template.vertex(i,1),1);
                if isempty(tmp)
                    Brain.vertex(i,1) = 0;
                else
                    Brain.vertex(i,1) = tmp;
                end
            end
            clearvars i tmp
            Brain.xyz = Brain.Template.xyz(Brain.Template2labels,:);
            Brain.labels = Brain.Template.labels(Brain.Template2labels,:);
            if isfield(Brain.Template,'regions')
                Brain.regions = Brain.Template.regions(1,Brain.Template2labels);
            end
            if isfield(Brain.Template,'regionsLR')
                Brain.regionsLR = Brain.Template.regionsLR(1,Brain.Template2labels);
            end
            disp(['   ' num2str(size(Brain.labels,1)) ' valid regions selected'])
            
            disp('   Calculate matrix parameters')
            Brain.matrixedges = lab_calc_edges(size(Brain.labels,1));
            
            progressbar('Match surface data to labels');
            disp('   Match surface data to labels')
            Vertex = zeros(size(Brain.vertices,1),1);
            for i = 1:size(Brain.Template2labels,1)
                Vertex(Brain.Template.vertex == Brain.Template2labels(i,1)) = i;
            end
            Vertex(Vertex==0) = i+1;
            tmp = Vertex(Brain.faces);
            Brain.mapsall = zeros(size(Brain.faces,1),size(Brain.Template2labels,1)+1);
            for j = 1:i+1
                Brain.mapsall(tmp(:,1) == j,j,1) = 1;
                Brain.mapsall(tmp(:,2) == j,j,2) = 1;
                Brain.mapsall(tmp(:,3) == j,j,3) = 1;
            end
            Brain.mapsall = mean(Brain.mapsall,3);
            clearvars i j Vertex tmp
            
            store = 1;
        end
    end
end

progressbar(1);

if store == 1 & isfield(Brain,'Brain_file')
    disp('   Save processed data to MNIbrain.mat')
    [~,Filepath] = lab_filename(Brain.Brain_file);
    save(fullfile(Filepath,'MNIbrain.mat'),'-struct','Brain');
    clearvars store Filepath
end

end

function labelsselect = load_labels(labelsselect,Hlabels)
   labels = get(Hlabels,'String');
   [Labels_file,Labels_filepath] = uigetfile('*.xls;*.xlsx','Select file with regions');
   if isnumeric(Labels_file)
       return
   end
   LabelsFile = lab_read_xls(fullfile(Labels_filepath,Labels_file));
   if isempty(LabelsFile)
       return
   end
   if size(LabelsFile,2) == 1
       LabelsFile = LabelsFile';
   end
   if verLessThan('matlab','8')
       [~,~,selection] = intersect(LabelsFile(1,:),labels);
       if isempty(selection)
           [~,~,selection] = intersect(LabelsFile(:,1),labels);
       end
   else
       [~,~,selection] = intersect(LabelsFile(1,:),labels,'stable');
       if isempty(selection)
           [~,~,selection] = intersect(LabelsFile(:,1),labels,'stable');
       end
   end
   if ~isempty(selection)
       labelsselect = labels(selection);
       labels = setdiff(labels,labelsselect);
       labels = cat(1,labelsselect(:),labels(:));
       set(Hlabels,'String',labels);
   end
end

function save_labels(labelsselect)
   [Labels_file,Labels_filepath] = uiputfile('SelectedRegions.xlsx','Select file to store regions');
   lab_write_xls(fullfile(Labels_filepath,Labels_file),labelsselect);
end