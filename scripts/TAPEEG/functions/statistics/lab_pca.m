% Calculate principal component analysis
%
% [Result,cfg] = lab_pca(data,header,cfg)
%
% data = array (subjects x measures)
% header = see lab_read_statistics
% cfg (optional)
%
% Written by F. Hatz 2014

function [Result,cfg] = lab_pca(data,header,cfg)

if ~exist('data','var')
    % read data
    [xlsout,~,~,~,cfg] = lab_read_statistics([],-1,0,0,1,0);
    if isempty(xlsout)
        Result = [];
        return
    end
    strlist = xlsout(2:end,1)';
    cfg.PCA.select = 1:length(strlist);
    cfg.PCA.method = 1;
    Prompt = {};
    Formats = {};
    Prompt(end+1,:) = {'Select measures','select'};
    Formats(end+1,1).type = 'list';
    Formats(end,1).style = 'listbox';
    Formats(end,1).items = strlist;
    Formats(end,1).limits = [0 4]; % multi-select
    Formats(end,1).size = [140 350];
    Formats(end+1,1).type = 'none';
    Prompt(end+1,:) = {'Find optimal number of factors','method'};
    Formats(end+1,1).type = 'list';
    Formats(end,1).style = 'popupmenu';
    Formats(end,1).items = {'Random','>80%','>90%','>95%','fixed'};
    [cfg.PCA,Cancelled] = inputsdlg(Prompt,'PCA',Formats,cfg.PCA);
    if Cancelled == 1 | isempty(cfg.PCA.select)
        return
    end
    data = cell2mat(xlsout(cfg.PCA.select+1,2:end))';
    header.vars = xlsout(cfg.PCA.select+1,1)';
    header.subjects = xlsout(1,2:end)';
    clearvars Prompt Formats
end
if ~exist('cfg','var')
    cfg = [];
end
if exist('header','var') & isfield(header,'vars')
    varnames = header.vars;
elseif exist('header','var') & isfield(header,'channels')
    if ~isfield(header,'numdatachannels')
        header.numdatachannels = size(data,1);
    end
    varnames = cellstr(header.channels(1:header.numdatachannels,:))';
    data = data(1:header.numdatachannels,:)';
else
    varnames = cellstr(num2str((1:size(data,2))'))';
end
if isfield(cfg,'file') & isfield(cfg,'path') & ~isempty(cfg.file) & ~isempty(cfg.path)
    [~,~,~,fileout] = lab_filename(cfg.file);
    warning off %#ok<WNOFF>
    mkdir(fullfile(cfg.path,'PCA'));
    warning off %#ok<WNOFF>
    fileout = fullfile(fullfile(cfg.path,'PCA'),fileout);
    cd(fullfile(cfg.path,'PCA'));
else
    fileout = [];
end

% calculate correlation matrix
[Result.C,Result.Cp] = corr(data,data);

% do pca
w = 1./var(data);
[wcoeff,score,latent,tsquared,explained] = pca(data,'VariableWeights',w);
coefforth = diag(sqrt(w))*wcoeff;
Result.PCA.wcoeff = wcoeff;
Result.PCA.coefforth = coefforth;
Result.PCA.score = score;
Result.PCA.latent = latent;
Result.PCA.tsquared = tsquared;
Result.PCA.explained = explained;
[~,index] = sort(tsquared,'descend'); % sort in descending order
Result.extreme = index(1);

if ~exist('cfg','var') | ~isfield(cfg,'PCA') | ~isfield(cfg.PCA,'method') | cfg.PCA.method == 1
    % do pca for randomized data (to find optimal number of components)
    rng('default');
    rng('shuffle');
    for j = 1:50
        for i = 1:size(data,2)
            tmp = randperm(size(data,1));
            datatmp(:,i) = data(tmp,i); %#ok<AGROW>
        end
        [~,~,~,~,Rexplained(:,j)] = pca(datatmp,'VariableWeights',w); %#ok<AGROW>
        clearvars tmp
    end
    for j = 1:50
        tmp(j) = find((explained-Rexplained(:,j)) > 0,1,'last');
    end
    nFactors = max(tmp);
    clearvars tmp
    Result.PCA.random_explained = Rexplained;
    Result.nFactors = nFactors;
elseif cfg.PCA.method == 2
    nFactors = 1;
    while sum(explained(1:nFactors)) < 80
        nFactors = nFactors + 1;
    end
    Result.nFactors = nFactors;
elseif cfg.PCA.method == 3
    nFactors = 1;
    while sum(explained(1:nFactors)) < 90
        nFactors = nFactors + 1;
    end
    Result.nFactors = nFactors;
elseif cfg.PCA.method == 4
    nFactors = 1;
    while sum(explained(1:nFactors)) < 95
        nFactors = nFactors + 1;
    end
    Result.nFactors = nFactors;
else
    answer = inputdlg('number of factors','Factors');
    nFactors = str2num(answer{1}); %#ok<ST2NM>
    Result.nFactors = nFactors;
end

% do factor analysis
warning off %#ok<WNOFF>
[Loadings,specificVar,T,stats,F] = factoran(data,nFactors,'rotate','promax');
warning on %#ok<WNON>
Result.FACTORS.coeff = Loadings;
Result.FACTORS.scores = F;
Result.FACTORS.specificvar = specificVar;
Result.FACTORS.T = T;
Result.FACTORS.stats = stats;

% Plot results
C = tril(Result.C) + triu(Result.Cp,1);
fig1 = figure('Color',[1 1 1],'NumberTitle','off','Name','Correlations','Menubar','none');
m1 = uimenu(fig1,'Label','File');
uimenu(m1,'Label','Save Picture','Callback',@(~,~)lab_print_figure);
uimenu(m1,'Label','Copy Picture','Callback','print -dmeta;');
uimenu(m1,'Label','Close','Callback','close;');
imagesc(Result.C);
cmap = flipud(gray(64));
set(gcf,'Colormap',cmap);
set(gca,'XTick',1:length(varnames),'XTickLabel',varnames,'YTick',1:length(varnames),'YTickLabel',varnames);
rotateXLabelsImage(gca,45);
axis normal
if ~isempty(fileout)
    lab_print_figure([fileout '_Correlations.pdf'],fig1);
end

fig2 = figure('Color',[1 1 1],'NumberTitle','off','Name','PC-Variance','Menubar','none');
m2 = uimenu(fig2,'Label','File');
uimenu(m2,'Label','Save Picture','Callback',@(~,~)lab_print_figure);
uimenu(m2,'Label','Copy Picture','Callback','print -dmeta;');
uimenu(m2,'Label','Close','Callback','close;');
pareto(explained);
xlabel('Principal Component');
ylabel('Variance Explained (%)');
if ~isempty(fileout)
    lab_print_figure([fileout '_PC_Variance.pdf'],fig2);
end

if nFactors > 1
    fig3 = figure('Color',[1 1 1],'NumberTitle','off','Name','Biplot-F1-F2','Menubar','none');
    m3 = uimenu(fig3,'Label','File');
    uimenu(m3,'Label','Save Picture','Callback',@(~,~)lab_print_figure);
    uimenu(m3,'Label','Copy Picture','Callback','print -dmeta;');
    uimenu(m3,'Label','Close','Callback','close;');
    biplot(Loadings(:,1:2),'scores',F(:,1:2),'varlabels',varnames);
    axis([-.26 0.6 -.51 .51]);
    if ~isempty(fileout)
        lab_print_figure([fileout '_Biplot_1_2.pdf'],fig3);
    end
end
if nFactors > 2
    fig4 = figure('Color',[1 1 1],'NumberTitle','off','Name','Biplot-F1-F3','Menubar','none');
    m4 = uimenu(fig4,'Label','File');
    uimenu(m4,'Label','Save Picture','Callback',@(~,~)lab_print_figure);
    uimenu(m4,'Label','Copy Picture','Callback','print -dmeta;');
    uimenu(m4,'Label','Close','Callback','close;');
    biplot(Loadings(:,[1 3]),'scores',F(:,[1 3]),'varlabels',varnames);
    axis([-.26 0.6 -.51 .51]);
    if ~isempty(fileout)
        lab_print_figure([fileout '_Biplot_1_3.pdf'],fig4);
    end
    
    fig5 = figure('Color',[1 1 1],'NumberTitle','off','Name','Biplot-F2-F3','Menubar','none');
    m5 = uimenu(fig5,'Label','File');
    uimenu(m5,'Label','Save Picture','Callback',@(~,~)lab_print_figure);
    uimenu(m5,'Label','Copy Picture','Callback','print -dmeta;');
    uimenu(m5,'Label','Close','Callback','close;');
    biplot(Loadings(:,2:3),'scores',F(:,2:3),'varlabels',varnames);
    axis([-.26 0.6 -.51 .51]);
    if ~isempty(fileout)
        lab_print_figure([fileout '_Biplot_2_3.pdf'],fig5);
    end
end

if exist('Rexplained','var')
    fig6 = figure('Color',[1 1 1],'NumberTitle','off','Name','Explained-RandomExplained','Menubar','none');
    m6 = uimenu(fig6,'Label','File');
    uimenu(m6,'Label','Save Picture','Callback',@(~,~)lab_print_figure);
    uimenu(m6,'Label','Copy Picture','Callback','print -dmeta;');
    uimenu(m6,'Label','Close','Callback','close;');
    T6 = plot(explained,'-k','LineWidth',1.5);
    hold on
    T6(end+1) = plot(mean(Rexplained,2),'-b','LineWidth',1.5);
    xlabel('Number of principal component');
    ylabel('Variance Explained (%)');
    legend(T6,{'explained','random'},'Location','NorthEast');
    if ~isempty(fileout)
        lab_print_figure([fileout '_explained_randomexpl.pdf'],fig6);
    end
end

fig7 = figure('Color',[1 1 1],'NumberTitle','off','Name','Coeff','Menubar','none');
m7 = uimenu(fig7,'Label','File');
uimenu(m7,'Label','Save Picture','Callback',@(~,~)lab_print_figure);
uimenu(m7,'Label','Copy Picture','Callback','print -dmeta;');
uimenu(m7,'Label','Close','Callback','close;');
T7 = bar(abs(coefforth)','stacked');
hold on
tmp = get(gca,'YTick');
tmp1 = explained / max(explained) * tmp(end);
T7(end+1) = plot(tmp1,'-k','LineWidth',1.5);
varnames7 = [varnames {'explained'}];
if exist('Rexplained','var')
    hold on
    tmp2 = min(Rexplained,[],2) / max(explained) * tmp(end);
    T7(end+1) = plot(tmp2,'-b','LineWidth',1.5);
    varnames7 = [varnames7 {'Rexplained'}];
end
legend(T7,varnames7,'Location','EastOutside');
pos = get(fig7,'Position');
pos(3) = pos(3)*1.4;
set(fig7,'Position',pos);
if ~isempty(fileout)
    lab_print_figure([fileout '_Coeff.pdf'],fig7);
end

fig8 = figure('Color',[1 1 1],'NumberTitle','off','Name','Factors','Menubar','none');
m8 = uimenu(fig8,'Label','File');
uimenu(m8,'Label','Save Picture','Callback',@(~,~)lab_print_figure);
uimenu(m8,'Label','Copy Picture','Callback','print -dmeta;');
uimenu(m8,'Label','Close','Callback','close;');
barh(Loadings');
legend(varnames,'Location','EastOutside');
for i = 1:size(Loadings,2)
    factorvar{1,i} = ['Factor ' num2str(i)]; %#ok<AGROW>
end
set(gca,'YDir','reverse','YTickLabel',factorvar);
pos = get(fig8,'Position');
pos(3) = pos(3)*1.4;
set(fig8,'Position',pos);
if ~isempty(fileout)
    lab_print_figure([fileout '_Factors.pdf'],fig8);
end

% write xls-output
lab_write_xls([fileout '_FactorWeights.xls'],Loadings);
lab_write_xls([fileout '_FactorScores.xls'],F);
lab_write_xls([fileout '_Explained.xls'],explained);
lab_write_xls([fileout '_WeightsCoefficients.xls'],coefforth);
lab_write_xls([fileout '_Correlations.xls'],Result.C);

% save matlab container
Result.vars = header.vars;
Result.subject = header.subjects;
save([fileout '_pca.mat'],'Result','cfg');

% write factors for analysis
if exist('xlsout','var') & exist('header','var') & isfield(header,'subjects')
    if strcmp(xlsout{end,1},'FileNr')
        for i = 1:length(header.subjects)
            tmp = strfind(header.subjects{i,1},'_');
            if ~isempty(tmp)
                Sprefix{i} = header.subjects{i,1}(1:tmp(1)-1); %#ok<AGROW>
                header.subjects{i,1} = header.subjects{i,1}(tmp(1)+1:end);
            end
            clearvars tmp
        end
        tmp = cell2mat(xlsout(end,2:end));
        [~,splitindex] = unique(tmp,'last');
        splitindex = [0 splitindex(:)'];
    else
        splitindex = [0 size(data,1)];
    end
    result = data * Loadings;
    for i = 1:size(result,2)
        factorvar{1,i} = ['F' num2str(i,'%02d')]; %#ok<AGROW>
    end
    resultvartmp = {''};
    xlsoutS = header.subjects(splitindex(1)+1:splitindex(2),:);
    for i = 1:length(splitindex)-1
        xlsoutS = [xlsoutS num2cell(result(splitindex(i)+1:splitindex(i+1),:))]; %#ok<AGROW>
        for j = 1:length(factorvar)
            if exist('Sprefix','var')
                resultvartmp{1,end+1} = [Sprefix{splitindex(i+1)} '_' factorvar{j}]; %#ok<AGROW>
            else
                resultvartmp{1,end+1} = factorvar{j}; %#ok<AGROW>
            end
        end
    end
    xlsoutS = cat(1,resultvartmp,xlsoutS)';
    fileout2 = [fileout '_ResultFactors'];
    if size(xlsoutS,2) < 255
        lab_write_xls([fileout2 '.xls'],xlsoutS);
    else
        lab_write_xls([fileout2 '.xlsx'],xlsoutS);
    end
end

end
