function [Brain,Mappings] = lab_brain_template2regions(Brain,Mappings,cfg)
    
if ~exist('cfg','var')
    cfg = [];
end
if ~exist('Mappings','var') | isempty(Mappings)
    Mappings = lab_load_mappings([],cfg,'MappingsIS.xls');
end
if isempty(Mappings) | ~isfield(Mappings,'mappingsChannels')
    return
end

if Mappings.mappingsChannels ~= length(Brain.labels)
    disp('Mappings not possible, mismatch in number of regions')
    return
end

Flag = false(1,length(Brain.labels)+1);
XYZ = zeros(length(Mappings.mappings),size(Brain.xyz,2));
MapsAll = zeros(size(Brain.mapsall,1),length(Mappings.mappings)+1);
Vertex = zeros(size(Brain.vertices,1),1);
for i = 1:length(Mappings.mappings)
    XYZ(i,:) = mean(Brain.xyz(Mappings.mappings{i},:),1); 
    MapsAll(:,i) = sum(Brain.mapsall(:,Mappings.mappings{i}),2);
    Flag(Mappings.mappings{i}) = true;
    for j = Mappings.mappings{i}
        Vertex(Brain.vertex==j) = i;
    end  
end
MapsAll(:,end) = sum(Brain.mapsall(:,Flag == false),2);
Brain.labels = Mappings.mappingstitle(:);
Brain.xyz = XYZ;
Brain.vertex = Vertex;
Brain.mapsall = MapsAll;
Brain.matrixedges = lab_calc_edges(size(Brain.labels,1));
Brain.Template2labels = [];
Brain.regions = [];
Brain.regionsLR = [];