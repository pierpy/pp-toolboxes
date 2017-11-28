% Calculate phase by searching for zero crossings, data between zero
% crossings is interpolated
%
% [phase,lag] = lab_zerocrossing(data)
%
% data    = matrix (channels,timepoints)
% phase   = matrix (channels,timepoints), data before first and after last
%                                         zero corss is omitted
% lag     = timepoint of first zero cross
%
% written by F. Hatz, Vumc Amsterdam 2013

function [phase,lag] = lab_zerocross(data,header,cfgfilt)

if exist('cfgfilt','var') & ~isempty(cfgfilt)
    data = lab_filter(data,header,cfgfilt,'novrb');
end

phase = zeros(size(data,1),size(data,2));
startE = 0;
stopE = size(data,2);
for i = 1:size(data,1)
    datatmp = data(i,:);
    datatmp = detrend(datatmp);
    phasetmp = zeros(1,size(data,2));
    x0 = datatmp(1);
    for j = 2:size(datatmp,2)
        if x0 <= 0 & datatmp(j) > 0
            if exist('x1','var')
                phasetmp(x1:j-1) = -pi:2*pi/(j-x1-1):pi;
            elseif j > startE
                startE = j;
            end
            x1 = j;
        end
        x0 = datatmp(j);
    end
    if exist('x1','var') & x1 <= stopE
        stopE = x1 - 1;
    end
    clearvars x0 j datatmp x1
    phase(i,:) = phasetmp;
end
clearvars i
phase = phase(:,startE:stopE);
lag = startE - 1;