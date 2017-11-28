% Helper file for lab_plot_IS
%
% Written by F. Hatz 2013

function lab_plot_distanceselection(include)

data = zeros(1,78);
data(include) = 1;
lab_plot_IS(data);