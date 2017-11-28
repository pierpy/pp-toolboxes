% function to do an inverse logit
%
% dataout = lab_alogit(datain)
%
% written by F. Hatz 2013

function dataout = lab_alogit(datain)

for s1 = 1:size(datain,1)
    for s2 = 1:size(datain,2)
        for s3 = 1:size(datain,3)
            dataout(s1,s2,s3) = exp(datain(s1,s2,s3)) / (1 + exp(datain(s1,s2,s3)));
        end
    end
end
clearvars s1 s2 s3
    
return