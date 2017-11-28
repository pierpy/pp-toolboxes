% Plot results (xls-file) as image
%
% lab_plot_xls2image(data,header)
%
%       data,header     output of lab_read_statistics
%
% Written by F. Hatz 2013

function lab_plot_xls2image(data,header)

if ~exist('data','var')
    [data,header,~,~,settings] = lab_read_statistics([],-1,0,1,0,1);
    if isempty(data)
        return
    end
end

if settings.clustervars > 1 & isfield(settings,'clustervars2') & settings.clustervars2 > 1
    data = reshape(data,size(data,1),settings.clustervars,settings.clustervars2,size(data,2)/(settings.clustervars*settings.clustervars2));
    data = permute(data,[1 4 2 3]);
    Measures = header.measures';
    Measures = regexprep(Measures,'_',' ');
    Vars = header.variables;
    Vars2 = header.variables2;
elseif settings.clustervars > 1
    data = reshape(data,size(data,1),settings.clustervars,size(data,2)/settings.clustervars);
    data = permute(data,[1 3 2 ]);
    Measures = header.measures';
    Measures = regexprep(Measures,'_',' ');
    Vars = header.variables;
    Vars2 = {};
else
    Measures = header.measures';
    Vars = {};
    Vars2 = {};
end
Subjects = header.subjects';
Subjects = regexprep(Subjects,'_',' ');

[settings,skipprocessing] = lab_set_plot_xls2image(settings,Subjects,Measures,Vars,Vars2);
if skipprocessing == 1
    return
end

data = data(settings.selectsubjects,settings.selectmeasures,settings.selectvars,settings.selectvars2);
data = permute(data,[settings.setsubjects settings.setmeasures settings.setvars settings.setvars2]);
Subjects = Subjects(1,settings.selectsubjects);
Measures = Measures(1,settings.selectmeasures);
if ~isempty(Vars)
    Vars = Vars(1,settings.selectvars);
end
if ~isempty(Vars2)
    Vars2 = Vars2(1,settings.selectvars2);
end
tmp = [{Subjects} {Measures} {Vars} {Vars2}];
tmp = tmp(1,[settings.setsubjects settings.setmeasures settings.setvars settings.setvars2]);
Subjects = tmp{1};
Measures = tmp{2};
Vars = tmp{3};
Vars2 = tmp{4};

tmp = strfind(settings.file,'_');
if ~isempty(tmp)
    Name = settings.file(1:tmp(end)-1);
else
    Name = settings.file;
end

Limits = [min(data(:)) max(data(:))];
if length(Measures) > 100
    Mskip = 4;
elseif length(Measures) > 75
    Mskip = 3;
elseif length(Measures) > 50
    Mskip = 2;
else
    Mskip = 1;
end
for i = 1:length(Measures)
    if mod(i,Mskip) ~= 0
        Measures{i} = ' ';
    end
end
if length(Subjects) > 100
    Skip = 4;
elseif length(Subjects) > 75
    Skip = 3;
elseif length(Subjects) > 50
    Skip = 2;
else
    Skip = 1;
end
for i = 1:length(Subjects)
    if mod(i,Skip) ~= 0
        Subjects{i} = ' ';
    end
end

fig1 = figure('Color',[1 1 1],'MenuBar','none','NumberTitle','off','Name',Name);
m1 = uimenu(fig1,'Label','File');
uimenu(m1,'Label','Save Picture','Callback',@(~,~)lab_print_figure);
uimenu(m1,'Label','Copy Picture','Callback','print -dmeta;');
uimenu(m1,'Label','Close','Callback','close;');
for N = 1:size(data,3)
    for M = 1:size(data,4)
        subplot(size(data,4),size(data,3),N+(M-1)*size(data,3));
        datatmp = data(:,:,N,M);
        Cmap = colormap('gray');
        if settings.flipcolormap == false
            Cmap = flipud(Cmap);
            colormap(Cmap);
        end
        imagesc(datatmp,Limits);
        set(gca,'XTick',1:length(Measures),'XTickLabel',Measures,'YTick',1:length(Subjects), ...
            'YTickLabel',Subjects,'TickLength',[0;0],'FontSize',settings.fontsize);
        hold on
        if ~isempty(Vars) & ~isempty(Vars2)
            title([Vars2{M} ' ' Vars{N}]);
        elseif ~isempty(Vars)
            title(Vars{N});
        end
        tmp = get(gca,'YLim');
        for i = 1.5:length(Measures)
            plot([i i],tmp,'Color',[0 0 0]);
        end
        tmp = get(gca,'XLim');
        for i = 1.5:length(Subjects)
            plot(tmp,[i i],'Color',[0 0 0]);
        end
        if length(unique(datatmp(:))) > 2 & settings.colorbar == true
            colorbar;
        end
        if length(Measures) > 10
            xticklabel_rotate90(1:Mskip:length(Measures),Measures(1:Mskip:length(Measures)));
        end
    end
end

end
