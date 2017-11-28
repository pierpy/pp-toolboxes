function lab_showmatrix(matrix,Mauto,Msave,savename)

Mnbr = size(matrix,3);
if ~exist('Msave','var')
    Msave = 0;
end
if ~exist('Mauto','var') & Mnbr > 1
    Mauto = 0.5;
end

f = figure;
colormap('gray');
cmap = colormap;
colormap(flipud(cmap));
for i = 1:Mnbr
    imagesc(matrix(:,:,i));
    if Msave == 1
        lab_print_figure([savename num2str(i) '.jpg'],f);
        clf;
    end
    if i<Mnbr
        pause(Mauto)
    end
end
if Msave == 1
    close;
end
