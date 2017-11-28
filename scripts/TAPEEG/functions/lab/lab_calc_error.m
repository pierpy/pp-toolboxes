function Result = lab_calc_error(data,header,cfg)

Result = [];
if ~exist('data','var')
    [data,header,~,~,cfg] = lab_read_statistics([],-1,0,1,0,1);
    if isempty(data)
        return
    end
end

Prompt = {'Select method','method'};
Formats.type = 'list';
Formats.style = 'popupmenu';
Formats.format = 'input';
Formats.items = {'Distance','Correlation'};
[cfg,Cancelled] = inputsdlg(Prompt,'Method',Formats,cfg);
if Cancelled == 1
    return
end

Mnum = size(data,2);
if strcmp(cfg.method,'Distance')
    fprintf('Calc error using distance')
    R = zeros(size(data,2),Mnum);
    for i = 1:Mnum
        R(:,i) = sum(abs(data - repmat(data(:,i),1,Mnum)));
        if mod(i,100) == 0
            fprintf('.')
        end
    end
    R(1:Mnum+1:end) = NaN;
    disp(':')
elseif strcmp(cfg.method,'Correlation')
    disp('Calc error using correlation')
    R = corr(data);
    R(1:Mnum+1:end) = NaN;
else
    return
end

if isfield(cfg,'clustervars') & cfg.clustervars > 1
    Rnum = floor(Mnum/cfg.clustervars);
    % R2 = zeros(Rnum,Rnum);
    % for i = 1:cfg.clustervars:Mnum
    %    for j = 1:cfg.clustervars:Mnum
    %        R2((i-1)/cfg.clustervars+1,(j-1)/cfg.clustervars+1) = nanmean(nanmean(R(i:i+cfg.clustervars-1,j:j+cfg.clustervars-1)));
    %    end
    % end
    R2 = zeros(Rnum,Rnum,cfg.clustervars);
    for k = 1:size(R2,3)
        for i = k:cfg.clustervars:Mnum
            for j = k:cfg.clustervars:Mnum
                R2((i-k)/cfg.clustervars+1,(j-k)/cfg.clustervars+1,k) = R(i,j);
            end
        end
    end
    R = mean(R2,3);
    Mnum = Rnum; %#ok<NASGU>
    clearvars R2 Rnum
end

if ~isfield(cfg,'clustervars2')
    cfg.clustervars2 = 1;
    if isempty(header,'variables2')
        header.variables2 = {''};
    end
end
Measures = {};
if isfield(header,'measures')
    for i = 1:length(header.measures)
        for j = 1:length(header.variables2)
            Measures{end+1,1} = [regexprep(header.measures{i},'_',' ') ' ' header.variables2{j}]; %#ok<AGROW>
        end
    end
else
    for i = 1:length(header.measures)
        for j = 1:length(header.variables2)
            Measures{end+1,1} = ['M' num2str(i) ' V' num2str(j)]; %#ok<AGROW>
        end
    end
end

Result.Error = R;
Result.Measures = Measures;

if cfg.clustervars2 > 1
    Rtmp = R;
    Rtmp(tril(true(size(R,2),size(R,2)),-1)) = NaN;
    NumM = length(header.measures);
    NumV = cfg.clustervars2;
    R2 = zeros(NumM,NumV);
    R2i = zeros(NumM,NumV);
    R3 = zeros(NumM,NumV);
    R3i = zeros(NumM,NumV);
    for i = 1:NumM
        for j = 1:NumV
            tmp1 = j:NumV:size(Rtmp,2);
            tmp2 = setdiff(1:size(Rtmp,2),tmp1);
            tmp1i = sum(~isnan(Rtmp((i-1)*NumV+j,tmp1)));
            tmp1 = nansum(Rtmp((i-1)*NumV+j,tmp1));
            R2(i,j) = R2(i,j) + tmp1;
            R2i(i,j) = R2i(i,j) + tmp1i;
            tmp2i = sum(~isnan(Rtmp((i-1)*NumV+j,tmp2)));
            tmp2 = nansum(Rtmp((i-1)*NumV+j,tmp2));
            R3(i,j) = R3(i,j) + tmp2;
            R3i(i,j) = R3i(i,j) + tmp2i;
        end
    end
    R2 = sum(R2,1) ./ sum(R2i,1);
    R3 = sum(R3,1) ./ sum(R3i,1);
    Result.IntraVar_Error = mean(R2);
    Result.InterVar_Error = mean(R3);
end

% create output_dir
cfg.path = fullfile(cfg.path,'CalcError');
warning off %#ok<WNOFF>
mkdir(cfg.path);
warning on %#ok<WNON>

% write xls
lab_write_xls(fullfile(cfg.path,[cfg.file '_Err.xlsx']),Result.Error);
if isfield(Result,'IntraVar_Error')
    xlsout = {'IntraVar_Error','InterVar_Error';Result.IntraVar_Error,Result.InterVar_Error};
    lab_write_xls(fullfile(cfg.path,[cfg.file '_ErrVar.xlsx']),xlsout);
end

% write Mat-Container
save(fullfile(cfg.path,[cfg.file '_Err.mat']),'Result');

% Plot figure
Fig = figure('NumberTitle','off','Name','Error''s','MenuBar','none','Color',[1 1 1]);
IM = imagesc(R);
set(gca,'XTick',1:length(Measures),'XTickLabel',Measures,'YTick',1:length(Measures),'YTickLabel',Measures,'TickLength',[0 0]);
xticklabel_rotate90(1:length(Measures),Measures);
colormap(flipud(gray))
colorbar
m1 = uimenu(Fig,'Label','File');
uimenu(m1,'Label','Save Picture','Callback',@(~,~)lab_print_figure);
uimenu(m1,'Label','Copy Picture','Callback','print -dmeta;');
uimenu(m1,'Label','Close','Callback','close;');
m2 = uimenu(Fig,'Label','Edit');
uimenu(m2,'Label','Set Min/Max','Callback',@(hObj,evd)setrange(IM));
lab_print_figure(fullfile(cfg.path,[cfg.file '.tif']),Fig);

end

function setrange(T)
    A = get(T,'Parent');
    CLim = get(A,'CLim');
    settings.minval = CLim(1);
    settings.maxval = CLim(2);
    Prompt = {'Min value','minval';'Max value','maxval'};
    Formats.type = 'edit';
    Formats.format = 'float';
    Formats.size = 30;
    Formats.limits = [-inf inf];
    Formats = [Formats Formats];
    [settings,Cancelled] = inputsdlg(Prompt,'Range',Formats,settings);
    if Cancelled ~= 1
        CLim = [settings.minval settings.maxval];
        set(A,'CLim',CLim);
    end
end
    