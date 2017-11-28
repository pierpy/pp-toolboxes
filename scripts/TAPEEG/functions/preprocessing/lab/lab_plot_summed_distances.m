function lab_plot_summed_distances
    
cfg.SEARCH.searchfolder = uigetdir('Select folder');
cfg.SEARCH.searchstring = {'*_FindDistance.xlsx'};

Files = lab_search_files(cfg);
if isempty(Files.Filelist)
    return
end

data = [];
for j = 1:length(Files.Filelist)
    datatmp = lab_read_xls(Files.Filelist{j},1);
    if size(datatmp,2) == 2
        try %#ok<TRYNC>
            if ischar(datatmp{1,1})
                datatmp = cell2mat(datatmp(2:end,:));
            else
                datatmp = cell2mat(datatmp);
            end
            data = cat(1,data,datatmp);
        end
    end
end

if isempty(data)
    return
end

DiffCluster = data(:,1)';
DiffNorm = data(:,2)';
x = DiffCluster(1:find(DiffCluster == max(DiffCluster),1,'first'));
pNorm = polyfit(DiffCluster,DiffNorm,5);
Fnorm = polyval(pNorm,x);
MinNorm = x(Fnorm==min(Fnorm));
Step = DiffCluster(2) - DiffCluster(1);
fig1 = figure('Color',[1 1 1],'NumberTitle','off','Name','Find optimal maximal distance','Menubar','none');
m1 = uimenu(fig1,'Label','File');
uimenu(m1,'Label','Save Picture','Callback',@(~,~)lab_print_figure);
uimenu(m1,'Label','Copy Picture','Callback','print -dmeta;');
uimenu(m1,'Label','Close','Callback','close;');
p1 = plot(DiffCluster,DiffNorm,'.');
set(p1,'Color',[0.5 0.5 0.5],'LineWidth',2)
hold on
p2 = plot(x,Fnorm,'-');
set(p2,'Color',[0.2 0.2 0.2],'LineWidth',2)
set(gca,'XLim',[x(1)-Step x(end)+Step]);
title(['Summed Error (minimal error: ' num2str(MinNorm) ')']);
xlabel('minimal electrode distance');
lab_print_figure(fullfile(cfg.SEARCH.searchfolder,'FindDistance_SummedError.pdf'),fig1);

end