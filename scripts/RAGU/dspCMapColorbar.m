function dspCMapColorbar(step,type)

cla
switch(type)
    case 'br'
        x = (-8:8) * step;
        cm = zeros(16,3);
        for i = 1:8
            cm(i,:)   = [((i-1) * 0.125) ((i-1) * 0.125)  1];
            cm(i+8,:) = [1 1-i*0.125 1-i*0.125];
        end
    otherwise
            error('Colormap not yet implemented');
end

axis([min(x) max(x) 0 1]);

for i = 1:size(cm,1)
    patch([x(i) x(i+1) x(i+1) x(i)],[0 0 1 1],cm(i,:));
end
axis on

set(gca,'YTickLabel',[]);        
        