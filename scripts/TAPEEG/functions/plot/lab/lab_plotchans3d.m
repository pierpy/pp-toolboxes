% Plot electrodes as 3d scatter
%
% lab_plotchans3d(loc_file,selection)
%
% Written by F. Hatz 2013

function  lab_plotchans3d(loc_file,selection)

if ischar(loc_file)
    loc_file = lab_read_locs(loc_file);
end
X = loc_file.x;
Y = loc_file.y;
Z = loc_file.z;
labels = loc_file.labels;

if ~exist('selection','var') | isempty(selection)
    selection = 1:length(X);
end

figure('color',[1 1 1]);

lim=1.05*max([X Y Z]);
eps=lim/20;
plot3(X,Y,Z,'ro')
hold on
plot3(X(selection),Y(selection),Z(selection),'b*')

plot3([0.08 0.12],[0 0],[0 0],'r','LineWidth',4) % nose
plot3([0 lim],[0 0],[0 0],'b--')                 % axis
plot3([0 0],[0 lim],[0 0],'g--')
plot3([0 0],[0 0],[0 lim],'r--')
plot3(0,0,0,'b+')
text(lim+eps,0,0,'X','HorizontalAlignment','center',...
	'VerticalAlignment','middle','Color',[0 0 0],...
	'FontSize',10);
text(0,lim+eps,0,'Y','HorizontalAlignment','center',...
	'VerticalAlignment','middle','Color',[0 0 0],...
	'FontSize',10);
text(0,0,lim+eps,'Z','HorizontalAlignment','center',...
	'VerticalAlignment','middle','Color',[0 0 0],...
	'FontSize',10);
for i = 1:length(selection)
    text(X(selection(i)),Y(selection(i)),Z(selection(i))+eps,labels(selection(i)), ...
        'HorizontalAlignment','center','VerticalAlignment','middle', ...
        'Color',[0 0 0],'FontSize',10);
end
set(gca,'visible','off','PlotBoxAspectRatio',[1 1 1],'xtick',[],'ytick',[],'DataAspectRatio',[1 1 1]);
rotate3d on

return
