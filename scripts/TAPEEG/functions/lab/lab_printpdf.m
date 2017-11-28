function lab_printpdf(Filename)

set(gcf, 'Units', 'centimeters')
set(gcf, 'PaperUnits','centimeters');
pos = get(gcf,'Position');
set(gcf, 'PaperSize', [pos(3) pos(4)]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 pos(3) pos(4)]);
print(gcf, '-dpdf', Filename)