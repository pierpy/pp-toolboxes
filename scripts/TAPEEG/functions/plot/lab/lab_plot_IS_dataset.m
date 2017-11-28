% Helper file for lab_plot_IS
%
% Written by F. Hatz 2013

function lab_plot_IS_dataset(data,subjects)

for i = 1:size(data,2)
    plot.data_file = [subjects{1,i} '_alpha1.tif'];
    plot.nosingle = 1;
    plot.backgroundcolor = [1 1 1];
    plot.textcolor = [0 0 0];
    plot.data_filepath = pwd;
    plot.Color1 = 'bluered';
    plot.MaxValue = 1;
    plot.MinValue = -1;
    datatmp = data(:,i)';
    datatmp = (datatmp - min(datatmp)) / (max(datatmp) - min(datatmp));
    lab_plot_IS(datatmp,plot);
end