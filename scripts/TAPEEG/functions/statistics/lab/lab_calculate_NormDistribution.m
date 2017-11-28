% Script to test data for Normal Distribution
%
% lab_calculate_NormDistribution(data,alpha,header)
%
%    data,header        see lab_read_statistics
%    alpha              default = 0.05
%
%
% Written by F. Hatz 2014

function lab_calculate_NormDistribution(data,alpha,header)

if ~exist('header','var') | ~isfield(header,'vars')
    for i = 1:size(data,2)
        header.vars{1,i} = ['Var' num2str(i)];
    end
    
end
if ~isfield(header,'file') | ~isfield(header,'path')
    fileout = 'TestNormalDistribution.xls';
    filepath = pwd;
else
    fileout = [header.file  'TestNormalDistribution.xls'];
    filepath = header.path;
end
if ~exist('alpha','var')
    alpha = 0.05;
end

warning off
mkdir(fullfile(filepath,'Q-Q-Plots'));
filepath = fullfile(filepath,'Q-Q-Plots');
warning on
fileout = fullfile(filepath,fileout);

for i = 1:size(data,2)
    datatmp = data(~isnan(data(:,i)),i);
    CDF = normcdf(datatmp, mean(datatmp), std(datatmp,1));
    [KolmogorovSmirnov(1,i),pKS(1,i)] = kstest(datatmp,[datatmp,CDF],alpha);
    [ShapiroWilk(1,i),pSW(1,i)] = swtest(datatmp,alpha);
    [LillieTest(1,i),pLillie(1,i)] = lillietest(datatmp,alpha);
    [AndersonDarling(1,i),pAD(1,i)] = AnDartest(datatmp,alpha);
    [JarqueBera(1,i),pJB(1,i)] = jbtest(datatmp,alpha);
    clearvars CDF datatmp
end

xlsout = cat(1,header.vars,num2cell([KolmogorovSmirnov;ShapiroWilk;LillieTest;AndersonDarling;JarqueBera]));
xlsout = cat(2,{'';'KolmogorovSmirnov';'ShapiroWilk';'Lilliefors';'AndersonDarling';'JarqueBera'},xlsout);
lab_write_xls(fileout,xlsout);
clearvars xlsout

fileout = [fileout(1:end-4) '_P.xls'];
xlsout = cat(1,header.vars,num2cell([pKS;pSW;pLillie;pAD;pJB]));
xlsout = cat(2,{'';'KolmogorovSmirnov';'ShapiroWilk';'Lilliefors';'AndersonDarling';'JarqueBera'},xlsout);
lab_write_xls(fileout,xlsout);
clearvars xlsout fileout

for i = 1:length(header.vars)
    vars{i} = regexprep(header.vars{i},{':',',',';','\','_'},' ');
end

for i = 1:size(data,2)
    f = figure('Visible','off','Color',[1 1 1]);
    qqplot(data(~isnan(data(:,i)),i));
    title(['Q-Q-Plot ' vars{i} ' vs. normal distribution'])
    lab_print_figure(fullfile(filepath,[header.vars{i} '.jpg']),f);
    close(f);
end

