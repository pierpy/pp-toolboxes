% Helper file for lab_plot_IS
%
% Written by F. Hatz 2013

function [edges1,edges2,edges3,edges,cfg] = lab_plot_create_edges(matrix,matrixedges,cfg)

% Create edges information
edges = (abs(matrix.matrix(matrixedges(:,3))) + abs(matrix.matrix(matrixedges(:,4)))) ./ 2;
mx = max(edges);
mn = min(edges);
if mx > mn
    edges = ((edges - mn) / (mx-mn))*0.5;
else
    edges = zeros(size(matrixedges,1),1);
end
clearvars mx mn

edges1 = (abs(matrix.matrix1(matrixedges(:,3))) + abs(matrix.matrix1(matrixedges(:,4)))) ./ 2;
mx = max(edges1);
mn = min(edges1);
if mx > mn
    edges1 = ((edges1 - mn) / (mx-mn))*0.5;
else
    edges1 = zeros(size(matrixedges,1),1);
end
clearvars mx mn

if ~isempty(cfg.PLOT.matrixmin2)
    edges2 = (abs(matrix.matrix2(matrixedges(:,3))) + abs(matrix.matrix2(matrixedges(:,4)))) ./ 2;
    mx = max(edges2);
    mn = min(edges2);
    if mx > mn
        edges2 = ((edges2 - mn) / (mx-mn))*0.5;
    else
        edges2 = zeros(size(matrixedges,1),1);
    end
    clearvars mx mn
else
    edges2 = [];
end

if ~isempty(cfg.PLOT.matrixmin3)
    edges3 = (abs(matrix.matrix3(matrixedges(:,3))) + abs(matrix.matrix3(matrixedges(:,4)))) ./ 2;
    mx = max(edges3);
    mn = min(edges3);
    if mx > mn
        edges3 = ((edges3 - mn) / (mx-mn))*0.5;
    else
        edges3 = zeros(size(matrixedges,1),1);
    end
    clearvars mx mn
else
    edges3 = [];
end