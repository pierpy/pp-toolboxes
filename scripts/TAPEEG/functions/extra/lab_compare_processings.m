% Script to compare results of two different processings using
% Pearson-correlations and ICC
%
% written by F. Hatz 2012

function lab_compare_processings

disp('Select Results auto')
[Rauto,header,~,~,cfg] = lab_read_statistics([],-1,0,1,1,1);
Rauto = Rauto';
if isfield(cfg,'numfiles') & cfg.numfiles > 1
    Ntrials = cfg.numfiles;
else
    answer = inputdlg('Number of trials');
    Ntrials = str2num(answer{1}); %#ok<ST2NM>
end
Nsubjects = size(Rauto,2) / Ntrials;
Nvars = size(Rauto,1);

disp('Select Results manual')
[Rmanual,~,~,~,cfgtmp] = lab_read_statistics([],0,0,1,1,1);
Rmanual = Rmanual';
if size(Rmanual,1) ~= size(Rauto,1) | size(Rmanual,2) ~= size(Rauto,2) | ...
        (isfield(cfgtmp,'numfiles') & cfgtmp.numfiles ~= Ntrials & cfgtmp.numfiles ~= 1)
    disp('Abort, mismatch of input')
    return
end
clearvars cfgtmp

settings.mode = 'A-1';
settings.bootstrap = 1000;
settings.alpha = 0.05;
tmp =  permute(reshape(Rauto,[Nvars Nsubjects Ntrials]),[1 3 2]);
for i = 1:Nsubjects
    S_ICCauto(i) = ICC(tmp(:,:,i),'A-1'); %#ok<AGROW>
end
tmp =  permute(reshape(Rauto,[Nvars Nsubjects Ntrials]),[2 3 1]);
for i = 1:Nvars
    V_ICCauto(i) = ICC(tmp(:,:,i),'A-1'); %#ok<AGROW>
    V_ICCauto2(:,i) = bootstrp(settings.bootstrap,@do_ICC,tmp(:,:,i),settings); %#ok<AGROW>
end
clearvars tmp

tmp =  permute(reshape(Rmanual,[Nvars Nsubjects Ntrials]),[1 3 2]);
for i = 1:Nsubjects
    S_ICCmanual(i) = ICC(tmp(:,:,i),'A-1'); %#ok<AGROW>
end
tmp =  permute(reshape(Rmanual,[Nvars Nsubjects Ntrials]),[2 3 1]);
for i = 1:Nvars
    V_ICCmanual(i) = ICC(tmp(:,:,i),'A-1'); %#ok<AGROW>
    V_ICCmanual2(:,i) = bootstrp(settings.bootstrap,@do_ICC,tmp(:,:,i),settings); %#ok<AGROW>
end
clearvars tmp

tmp = corr(Rauto,Rmanual);
CORRsubjects = diag(tmp);

tmp = corr(Rauto',Rmanual');
CORRvars = diag(tmp);

for i = 1:Ntrials
    Start = (i-1)*Nsubjects + 1;
    Stop = i * Nsubjects;
    tmp = corr(Rauto(:,Start:Stop)',Rmanual(:,Start:Stop)');
    CORRvarsS(:,i) = diag(tmp); %#ok<AGROW>
end

[File,Path] = uiputfile('*.xlsx','Select file to store');
[~,~,~,File] = lab_filename(File);
Fileout = fullfile(Path,File);

for i = 1:Nsubjects
    tmp = strfind(header.subjects{i,1},'_');
    if ~isempty(tmp)
        subjects{i,1} = header.subjects{i,1}(1:tmp(end)-1); %#ok<AGROW>
    end
end
xlsout = [subjects num2cell(S_ICCmanual') num2cell(S_ICCauto')];
xlsout = cat(1,{'','ICC_manual','ICCauto'},xlsout);
lab_write_xls([Fileout '_ICCsubject.xlsx'],xlsout);

xlsout = [header.vars' num2cell(V_ICCmanual') num2cell(V_ICCauto')];
xlsout = cat(1,{'','ICC_manual','ICCauto'},xlsout);
lab_write_xls([Fileout '_ICCmeasure.xlsx'],xlsout);

xlsout = cat(1,header.vars,num2cell(V_ICCauto2));
lab_write_xls([Fileout '_ICCmeasure_auto_bootstrap.xlsx'],xlsout);

xlsout = cat(1,header.vars,num2cell(V_ICCmanual2));
lab_write_xls([Fileout '_ICCmeasure_manual_bootstrap.xlsx'],xlsout);

xlsout = [header.subjects num2cell(CORRsubjects)];
lab_write_xls([Fileout '_CORRsubject.xlsx'],xlsout);

xlsout = [header.vars' num2cell(CORRvars)];
lab_write_xls([Fileout '_CORRmeasure.xlsx'],xlsout);

xlsout = {'','Mean','Std'};
if isfield(cfg,'clustervars')
    clusters = cfg.clustervars;
else
    clusters = 1;
end
numvals = size(CORRvarsS,1) / clusters;
for i = 1:size(CORRvarsS,2)
    for j = 1:numvals
        Start = (j-1)*clusters+1;
        Stop = j*clusters;
        tmp = strfind(header.vars{1,Start},'_');
        if ~isempty(tmp)
            xlsout{end+1,1} = ['Trial' num2str(i) '_' header.vars{1,Start}(1:tmp(1)-1)]; %#ok<AGROW>
        else
            xlsout{end+1,1} = ['Trial' num2str(i) '_' header.vars{1,Start}]; %#ok<AGROW>
        end
        xlsout{end,2} = mean(CORRvarsS(Start:Stop,i));
        xlsout{end,3} = std(CORRvarsS(Start:Stop,i));
    end
end
lab_write_xls([Fileout '_CORRmeasureS.xlsx'],xlsout);

tmp = reshape(Rauto,[Nvars Nsubjects Ntrials]);
numel = 1;
for i = 1:Ntrials
    for j = 1:Ntrials
        if i~=j
            tmp2(:,:,numel) = abs(tmp(:,:,i) - tmp(:,:,j)); %#ok<AGROW>
            numel = numel+1;
        end
    end
end
ImageAuto = mean(tmp,3);
clearvars tmp numel i j tmp2

tmp = reshape(Rmanual,[Nvars Nsubjects Ntrials]);
numel = 1;
for i = 1:Ntrials
    for j = 1:Ntrials
        if i~=j
            tmp2(:,:,numel) = abs(tmp(:,:,i) - tmp(:,:,j)); %#ok<AGROW>
            numel = numel+1;
        end
    end
end
ImageManual = mean(tmp,3);
clearvars tmp numel i j tmp2

tmp = max([ImageAuto ImageManual],[],2);
ImageAuto2 = ImageAuto ./ repmat(tmp,1,size(ImageAuto,2));
ImageManual2 = ImageManual ./ repmat(tmp,1,size(ImageAuto,2));

MaxVal = max(max(ImageAuto(:)),max(ImageManual(:)));
Fontsize = 300/Nvars;
if Fontsize > 10
    Fontsize = 10;
end

fig1 = figure('Color',[1 1 1],'NumberTitle','off','Name','Kmeans-Clustering -- Normalized variance','Menubar','none');
m1 = uimenu(fig1,'Label','File');
uimenu(m1,'Label','Save Picture','Callback',@(~,~)lab_print_figure);
uimenu(m1,'Label','Copy Picture','Callback','print -dmeta;');
uimenu(m1,'Label','Close','Callback','close;');
imagesc(ImageAuto,[0 MaxVal]);
title('Mean Differences Auto');
xlabel('Subjects');
set(gca,'YTick',1:Nvars,'YTickLabel',header.vars,'Fontsize',Fontsize);
colorbar
lab_print_figure([Fileout '_DiffAuto.tiff'],fig1);

fig2 = figure('Color',[1 1 1],'NumberTitle','off','Name','Kmeans-Clustering -- Normalized variance','Menubar','none');
m1 = uimenu(fig1,'Label','File');
uimenu(m1,'Label','Save Picture','Callback',@(~,~)lab_print_figure);
uimenu(m1,'Label','Copy Picture','Callback','print -dmeta;');
uimenu(m1,'Label','Close','Callback','close;');
imagesc(ImageManual,[0 MaxVal]);
title('Mean Differences Manual');
xlabel('Subjects');
set(gca,'YTick',1:Nvars,'YTickLabel',header.vars,'Fontsize',Fontsize);
colorbar
lab_print_figure([Fileout '_DiffManual.tiff'],fig2);

fig3 = figure('Color',[1 1 1],'NumberTitle','off','Name','Kmeans-Clustering -- Normalized variance','Menubar','none');
m1 = uimenu(fig1,'Label','File');
uimenu(m1,'Label','Save Picture','Callback',@(~,~)lab_print_figure);
uimenu(m1,'Label','Copy Picture','Callback','print -dmeta;');
uimenu(m1,'Label','Close','Callback','close;');
imagesc(ImageAuto2,[0 1]);
title('Mean Differences Auto (normalized)');
xlabel('Subjects');
set(gca,'YTick',1:Nvars,'YTickLabel',header.vars,'Fontsize',Fontsize);
lab_print_figure([Fileout '_DiffAutoNorm.tiff'],fig3);

fig4 = figure('Color',[1 1 1],'NumberTitle','off','Name','Kmeans-Clustering -- Normalized variance','Menubar','none');
m1 = uimenu(fig1,'Label','File');
uimenu(m1,'Label','Save Picture','Callback',@(~,~)lab_print_figure);
uimenu(m1,'Label','Copy Picture','Callback','print -dmeta;');
uimenu(m1,'Label','Close','Callback','close;');
imagesc(ImageManual2,[0 1]);
title('Mean Differences Manual (normalized)');
xlabel('Subjects');
set(gca,'YTick',1:Nvars,'YTickLabel',header.vars,'Fontsize',Fontsize);
lab_print_figure([Fileout '_DiffManualNorm.tiff'],fig4);

fig5 = figure('Color',[1 1 1],'NumberTitle','off','Name','Kmeans-Clustering -- Normalized variance','Menubar','none');
m1 = uimenu(fig1,'Label','File');
uimenu(m1,'Label','Save Picture','Callback',@(~,~)lab_print_figure);
uimenu(m1,'Label','Copy Picture','Callback','print -dmeta;');
uimenu(m1,'Label','Close','Callback','close;');
boxplot(V_ICCauto2);
title('Boxplot ICC''s Auto');
lab_print_figure([Fileout '_Boxplot_ICC_Auto.tiff'],fig5);

fig6 = figure('Color',[1 1 1],'NumberTitle','off','Name','Kmeans-Clustering -- Normalized variance','Menubar','none');
m1 = uimenu(fig1,'Label','File');
uimenu(m1,'Label','Save Picture','Callback',@(~,~)lab_print_figure);
uimenu(m1,'Label','Copy Picture','Callback','print -dmeta;');
uimenu(m1,'Label','Close','Callback','close;');
boxplot(V_ICCmanual2);
title('Boxplot ICC''s Manual');
lab_print_figure([Fileout '_Boxplot_ICC_Manual.tiff'],fig6);

end

function Result = do_ICC(data,settings)
   Result = ICC(data,settings.mode,settings.alpha);
end
