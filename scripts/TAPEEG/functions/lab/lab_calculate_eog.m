function settings = lab_calculate_eog(data,header,settings)

if ~isfield(settings,'eog') | isempty(settings.eog)
    settings.zvalueeye = [];
    disp('    no eog could be found/calculated')
    return
end

if isfield(header,'eog') & ~isempty(header.eog)
    EOG = header.eog(:);
else
    EOG = settings.eog;
end
settings.blinkchans = [];
for n = 1:size(EOG,1)
    if EOG(n,1) > 0
        tmp(1) = EOG(n,1);
        if size(EOG,2) > 1 & EOG(n,2) > 0
            tmp(2) = EOG(n,2);
        end
        if size(EOG,2) > 2 & EOG(n,3) > 0 & isfield(header,'bad') & ...
                isfield(header.bad,'All') & find(header.bad.All==EOG(n,1))
            tmp(1) = EOG(n,3);
        end
        if size(EOG,2) > 3 & EOG(n,4) > 0 & isfield(header,'bad') & ...
                isfield(header.bad,'All') & find(header.bad.All==EOG(n,2))
            tmp(2) = EOG(n,4);
        end
        if length(tmp) == 1 & tmp(1) <= size(data,1)
            settings.blinkchans = cat(1,settings.blinkchans,data(tmp(1),:));
        elseif length(tmp) == 2 & tmp(1) <= size(data,1) & tmp(2) <= size(data,1)
            settings.blinkchans = cat(1,settings.blinkchans,data(tmp(1),:) - data(tmp(2),:));
        end
    end
end