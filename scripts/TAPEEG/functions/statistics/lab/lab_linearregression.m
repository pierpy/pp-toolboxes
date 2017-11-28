% Script to correct for influence of one or multiple factors on data.
%
% [dataR,stats] = lab_linearregression(data,factors,result)
%
% Output is a matrix with residuals
%
% data    = matrix(subjects x variables)
% factors = matrix(subjects x variables)
%
% Written by F. Hatz 2012 Neurology University Hospital Basel

function [dataR,stats] = lab_linearregression(data,factors,result)
if islogical(data)
    distribution = 'binominal';
else
    distribution = 'normal';
end
disp(['Generalized linear model regression to correct for factors (distribution: ' distribution ')'])

for Nvar = 1:size(data,2)
    [~,~,statstmp] = glmfit([factors result],data(:,Nvar),distribution);
    dataR(:,Nvar) = statstmp.resid + (result .* statstmp.beta(end)); %#ok<AGROW>
    VarsNr = size(statstmp.beta,1);
    statsR(1,Nvar) = statstmp.dfe; %#ok<AGROW>
    statsR(2,Nvar) = statstmp.sfit; %#ok<AGROW>
    statsR(3:2+VarsNr,Nvar) = statstmp.beta; %#ok<AGROW>
    statsR(3+VarsNr:2+2*VarsNr,Nvar) = statstmp.t; %#ok<AGROW>
    statsR(3+2*VarsNr:2+3*VarsNr,Nvar) = statstmp.p; %#ok<AGROW>
end

stats = cell(3 * VarsNr + 2,size(data,2)+1);
stats{1,1} = 'dfe';
stats{2,1} = 'sfit';
for i = 1:VarsNr
    stats{(2+i),1} = ['beta_' num2str(i)];
    stats{(2+i+VarsNr),1} = ['T_' num2str(i)];
    stats{(2+i+VarsNr*2),1} = ['pvalue_' num2str(i)];
end
stats(:,2:end) = num2cell(statsR);
