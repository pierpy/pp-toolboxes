% Result = lab_pca2clustering(data)
%
% rows are patients / trials
% columns are variables per patient/trial

function Result = lab_matrix2clustering(matrix,noplots)

disp('   Do PCA to reduce number of links')

if ~exist('noplots','var')
    noplots = false;
end

nchans = size(matrix,1);
matrix = reshape(matrix,nchans^2,size(matrix,3))';
index = lab_get_lowertriangle(nchans,1);
matrix = matrix(:,index);

% calculate correlation matrix
Result.C = corr(matrix,matrix);

% do pca
w = 1./var(matrix);
[wcoeff,score,latent,tsquared,explained] = pca(matrix,'VariableWeights',w);
coefforth = diag(sqrt(w))*wcoeff;
Result.PCA.wcoeff = wcoeff;
Result.PCA.coefforth = coefforth;
Result.PCA.score = score;
Result.PCA.latent = latent;
Result.PCA.tsquared = tsquared;
Result.PCA.explained = explained;
[~,index] = sort(tsquared,'descend'); % sort in descending order
Result.extreme = index(1);

rng('default');
rng('shuffle');
for j = 1:50
    for i = 1:size(matrix,1)
        tmp = randperm(size(matrix,2));
        matrixtmp(i,:) = matrix(i,tmp);
    end
    [~,~,~,~,Rexplained(:,j)] = pca(matrixtmp,'VariableWeights',w);
    clearvars tmp
end
for j = 1:50
    tmp(j) = find((explained-Rexplained(:,j)) > 0,1,'last');
end
nFactors = max(tmp);
clearvars tmp
Result.PCA.random_explained = Rexplained;
Result.nFactors = nFactors;
disp(['    Number of components set to ' num2str(nFactors)])

if noplots == false
    varnames = cellstr(num2str((1:nchans)'));
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
    
    fig2 = figure('Color',[1 1 1],'NumberTitle','off','Name','PC-Variance','Menubar','none');
    m2 = uimenu(fig2,'Label','File');
    uimenu(m2,'Label','Save Picture','Callback',@(~,~)lab_print_figure);
    uimenu(m2,'Label','Copy Picture','Callback','print -dmeta;');
    uimenu(m2,'Label','Close','Callback','close;');
    pareto(explained);
    xlabel('Principal Component');
    ylabel('Variance Explained (%)');
    
    fig3 = figure('Color',[1 1 1],'NumberTitle','off','Name','Explained-RandomExplained','Menubar','none');
    m3 = uimenu(fig3,'Label','File');
    uimenu(m3,'Label','Save Picture','Callback',@(~,~)lab_print_figure);
    uimenu(m3,'Label','Copy Picture','Callback','print -dmeta;');
    uimenu(m3,'Label','Close','Callback','close;');
    T3 = plot(explained,'-k','LineWidth',1.5);
    hold on
    T3(end+1) = plot(mean(Rexplained,2),'-b','LineWidth',1.5);
    xlabel('Number of principal component');
    ylabel('Variance Explained (%)');
    legend(T3,{'explained','random'},'Location','NorthEast');
end

end