% Read xls files with channel mappings info
%
% Mappings = lab_read_mappings(cfg,ISflag)
%
% cfg    = structure with config (optional)
% ISflag = as default the file 'MappingsIS.xls' is loaded
%          (in place of 'Mappings.xls')
%
% written by F. Hatz 2012

function Mappings = lab_read_mappings(mappings_file)

if ~exist('mappings_file','var')
    [mappings_file,mappings_filepath] = uigetfile('*.xls;*.xlsx','Select Mappings-file');
else
    [mappings_file,mappings_filepath] = lab_filename(mappings_file);
end

if ~exist(fullfile(mappings_filepath,mappings_file),'file')
    Mappings = [];
    return
end

try
    disp(['   Read (' mappings_file ')'])
    if ispc
        [~,~,mappingstmp] = xlsread(fullfile(mappings_filepath,mappings_file),1);
    else
        [~,~,mappingstmp] = xlsread(fullfile(mappings_filepath,mappings_file),1,'','basic');
    end
catch %#ok<CTCH>
    Mappings = [];
    return
end

if size(mappingstmp,2) == 2
    mappingstitle = mappingstmp(:,1);
    for i = 1:size(mappingstmp,1)
        if isnumeric(mappingstmp{i,2})
            mtmp1{i} = mappingstmp{i,2}; %#ok<AGROW>
        else
            mtmp1{i} = str2num(mappingstmp{i,2}); %#ok<ST2NM,AGROW>
        end
        if ~isempty(mtmp1{i}) & ~isnan(mtmp1{i})
            mtmp2(mtmp1{i}) = i; %#ok<AGROW>
        end
    end
elseif strcmp(mappingstmp{1,1},'Region')
    mappingstitle = mappingstmp(2:end,1);
    mappingstmp = cell2mat(mappingstmp(2:end,2:end));
    for i = 1:size(mappingstmp,1)
        mtmp1{i} = find(mappingstmp(i,:)>0); %#ok<AGROW>
        mtmp2(mappingstmp(i,:)>0) = i; %#ok<AGROW>
    end
end
if ~exist('mtmp1','var') | ~exist('mtmp2','var')
    Mappings = [];
    return
end
Mappings.mappingsChannelsFile = length(mtmp2);
Mappings.mappingsChannels = length(mtmp2);
Mappings.mappings = mtmp1;
Mappings.mappingstitle = mappingstitle;
clearvars mappingstitle mappingstmp i
if strcmp(Mappings.mappingstitle{end,1},'no region') | isempty(Mappings.mappings{1,end})
    Mappings.mappingstitle = Mappings.mappingstitle(1:end-1,1);
    Mappings.mappings = Mappings.mappings(1,1:end-1);
end
for i = 1:size(Mappings.mappingstitle,1)
    tmp = Mappings.mappingstitle{i,1};
    tmp=textscan(tmp,'%s');
    tmp = tmp{1,1};
    title = '';
    for j = 1:size(tmp,1);
        if size(tmp,1) == 1 & length(tmp{1,1}) > 1
            title = [upper(tmp{1,1}(1)) lower(tmp{1,1}(2))];
        elseif j < 4
            title = [title upper(tmp{j,1}(1))];
        end
    end
    Mappings.mappingstitleS{i,1} = title;
end
