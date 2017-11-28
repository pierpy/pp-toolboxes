function dspPMapColorbar(LevelList)
cla

nNegContours = sum(LevelList < 0);
nPosContours = sum(LevelList > 0) ;               

axis([0 (nNegContours+nPosContours) 0 1]);


for i = 1:nNegContours
    colc = (i-1) / nNegContours;
    col = [colc colc 1];             
    patch([i-1 i i i-1],[0 0 1 1],col);
end

for j = 1:nPosContours
    colc = j / nPosContours;
    col = [1 1-colc 1-colc];
    i = j + nNegContours;
    patch([i-1 i i i-1],[0 0 1 1],col);
end

axis on

set(gca,'YTickLabel',[],'XTick',0:(nNegContours+nPosContours),'XTickLabel',abs(LevelList));
xticklabel_rotate([],90,[],'Fontsize',6);        