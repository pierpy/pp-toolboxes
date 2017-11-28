function lab_calc_matrixcorr

[data,header,~,~,cfg] = lab_read_statistics([],0,0,1,0,1);
if ~isfield(cfg,'clustervars2')
    cfg.clustervars2 = 1;
end
if ~isfield(cfg,'clustervars')
    cfg.clustervars = 1;
end
cfg.nummeasures = (size(data,2)/(cfg.clustervars*cfg.clustervars2));
dat = reshape(data,[size(data,1) cfg.clustervars cfg.clustervars2 cfg.nummeasures]);

R1 = zeros(size(dat,3),size(dat,4));
P1 = zeros(size(dat,3),size(dat,4));
for i = 1:size(dat,3)
    for j = 1:size(dat,4)
        [tmp,ptmp] = corr(dat(:,:,i,j));
        R1(i,j) = mean(tmp(tril(true(size(tmp,1),size(tmp,2)),-1)));
        P1(i,j) = mean(ptmp(tril(true(size(tmp,1),size(tmp,2)),-1)));
        clearvars tmp ptmp
    end
end

dat2 = reshape(dat,[size(dat,1),size(dat,2),size(dat,3)*size(dat,4)]);
R2 = zeros(size(dat2,3),size(dat2,3));
P2 = zeros(size(dat2,3),size(dat2,3));
for i = 1:size(dat2,3)
    for j = 1:size(dat2,3)
        tmp = zeros(1,size(dat2,1));
        ptmp = zeros(1,size(dat2,1));
        for m = 1:size(dat2,1)
            [tmp2,ptmp2] = corr(permute(dat2(1,:,[i j]),[2 3 1]));
            tmp(m) = tmp2(2,1);
            ptmp(m) = ptmp2(2,1);
        end
        R2(i,j) = mean(tmp);
        P2(i,j) = mean(ptmp);
        clearvars tmp ptmp tmp2 ptmp2
    end
end

R3 = reshape(R2,[cfg.clustervars2 cfg.nummeasures cfg.clustervars2 cfg.nummeasures]);
R3 = permute(R3,[2 1 4 3]);
R3 = reshape(R3,[size(R2,1) size(R2,2)]);

P3 = reshape(P2,[cfg.clustervars2 cfg.nummeasures cfg.clustervars2 cfg.nummeasures]);
P3 = permute(P3,[2 1 4 3]);
P3 = reshape(P3,[size(P2,1) size(P2,2)]);

R1(1:size(R1,1)+1:end) = NaN;
R2(1:size(R2,1)+1:end) = NaN;
R3(1:size(R3,1)+1:end) = NaN;

P1(P1>0.05) = 0.05;
P2(P2>0.05) = 0.05;
P3(P3>0.05) = 0.05;

P1(1:size(P1,1)+1:end) = NaN;
P2(1:size(P2,1)+1:end) = NaN;
P3(1:size(P3,1)+1:end) = NaN;

fig1 = figure('Color',[1 1 1],'NumberTitle','off','Name','Pearson-Correlation','Menubar','none');
m1 = uimenu(fig1,'Label','File');
uimenu(m1,'Label','Save Picture','Callback',@(~,~)lab_print_figure);
uimenu(m1,'Label','Copy Picture','Callback','print -dmeta;');
uimenu(m1,'Label','Close','Callback','close;');
imagesc(R1);
title('Mean correlations');
set(gca,'YTick',1:cfg.clustervars2,'YTickLabel',header.variables2);
set(gca,'XTick',1:cfg.nummeasures,'XTickLabel',header.measures);
colorbar;
colormap(flipud(gray));

fig1p = figure('Color',[1 1 1],'NumberTitle','off','Name','P-Value','Menubar','none');
m1p = uimenu(fig1p,'Label','File');
uimenu(m1p,'Label','Save Picture','Callback',@(~,~)lab_print_figure);
uimenu(m1p,'Label','Copy Picture','Callback','print -dmeta;');
uimenu(m1p,'Label','Close','Callback','close;');
imagesc(P1);
title('Mean P-Values');
set(gca,'YTick',1:cfg.clustervars2,'YTickLabel',header.variables2);
set(gca,'XTick',1:cfg.nummeasures,'XTickLabel',header.measures);
colorbar;
colormap(gray)

allmeasures = {};
for i = 1:cfg.nummeasures
    for j = 1:cfg.clustervars2
        allmeasures{end+1} = [header.measures{i} ' ' header.variables2{j}]; %#ok<AGROW>
    end
end
fig2 = figure('Color',[1 1 1],'NumberTitle','off','Name','Pearson-Correlation','Menubar','none');
m2 = uimenu(fig2,'Label','File');
uimenu(m2,'Label','Save Picture','Callback',@(~,~)lab_print_figure);
uimenu(m2,'Label','Copy Picture','Callback','print -dmeta;');
uimenu(m2,'Label','Close','Callback','close;');
imagesc(R2);
title('Mean correlations');
set(gca,'YTick',1:length(allmeasures),'YTickLabel',allmeasures);
set(gca,'XTick',1:length(allmeasures),'XTickLabel',allmeasures);
xticklabel_rotate90(1:length(allmeasures),allmeasures);
colorbar;
colormap(flipud(gray));

fig2p = figure('Color',[1 1 1],'NumberTitle','off','Name','P-Value','Menubar','none');
m2p = uimenu(fig2p,'Label','File');
uimenu(m2p,'Label','Save Picture','Callback',@(~,~)lab_print_figure);
uimenu(m2p,'Label','Copy Picture','Callback','print -dmeta;');
uimenu(m2p,'Label','Close','Callback','close;');
imagesc(P2);
title('Mean P-Values');
set(gca,'YTick',1:length(allmeasures),'YTickLabel',allmeasures);
set(gca,'XTick',1:length(allmeasures),'XTickLabel',allmeasures);
xticklabel_rotate90(1:length(allmeasures),allmeasures);
colorbar;
colormap(gray)

allmeasures2 = {};
for i = 1:cfg.clustervars2
    for j = 1:cfg.nummeasures
        allmeasures2{end+1} = [header.variables2{i} ' ' header.measures{j}]; %#ok<AGROW>
    end
end
fig3 = figure('Color',[1 1 1],'NumberTitle','off','Name','Pearson-Correlation','Menubar','none');
m3 = uimenu(fig3,'Label','File');
uimenu(m3,'Label','Save Picture','Callback',@(~,~)lab_print_figure);
uimenu(m3,'Label','Copy Picture','Callback','print -dmeta;');
uimenu(m3,'Label','Close','Callback','close;');
imagesc(R3);
title('Mean correlations');
set(gca,'YTick',1:length(allmeasures2),'YTickLabel',allmeasures2);
set(gca,'XTick',1:length(allmeasures2),'XTickLabel',allmeasures2);
xticklabel_rotate90(1:length(allmeasures2),allmeasures2);
colorbar;
colormap(flipud(gray));

fig3p = figure('Color',[1 1 1],'NumberTitle','off','Name','P-Value','Menubar','none');
m3p = uimenu(fig3p,'Label','File');
uimenu(m3p,'Label','Save Picture','Callback',@(~,~)lab_print_figure);
uimenu(m3p,'Label','Copy Picture','Callback','print -dmeta;');
uimenu(m3p,'Label','Close','Callback','close;');
imagesc(P3);
title('Mean P-Values');
set(gca,'YTick',1:length(allmeasures2),'YTickLabel',allmeasures2);
set(gca,'XTick',1:length(allmeasures2),'XTickLabel',allmeasures2);
xticklabel_rotate90(1:length(allmeasures2),allmeasures2);
colorbar;
colormap(gray)

lab_print_figure(fullfile(cfg.path,[cfg.file '_Subjects.emf']),fig1);
lab_print_figure(fullfile(cfg.path,[cfg.file '_SubjectsP.emf']),fig1p);
lab_print_figure(fullfile(cfg.path,[cfg.file '_Single1.emf']),fig2);
lab_print_figure(fullfile(cfg.path,[cfg.file '_Single1P.emf']),fig2p);
lab_print_figure(fullfile(cfg.path,[cfg.file '_Single2.emf']),fig3);
lab_print_figure(fullfile(cfg.path,[cfg.file '_Single2P.emf']),fig3p);