% Helper file for lab_test_quality
%
% written by F. Hatz Vumc 2013

function R = lab_evaluate_quality(R,settings)

if isempty(settings) | isempty(R)
    R.good = false;
    return
end

R.bchans = [];
tmp = 0;
if isfield(settings,'cogfreq') & ~isempty(settings.cogfreq) & ...
        isfield(R,'cogfreq') & R.cogfreq > settings.cogfreq
    tmp = tmp + 1;
end
if isfield(settings,'peakfreq') & ~isempty(settings.peakfreq) & ...
        isfield(R,'peakfreq') & R.peakfreq > settings.peakfreq
    tmp = tmp + 1;
end
if isfield(settings,'peakamp') & ~isempty(settings.peakamp) & ...
        isfield(R,'peakamp') & R.peakamp > settings.peakamp
    tmp = tmp + 1;
end
if isfield(settings,'peak2min') & ~isempty(settings.peak2min) & ...
        isfield(R,'peak2min') & R.peak2min > settings.peak2min
    tmp = tmp + 1;
end
if isfield(settings,'peakratio') & ~isempty(settings.peakratio) & ...
        isfield(R,'peakratio') & R.peakratio > settings.peakratio
    tmp = tmp + 1;
end
if isfield(settings,'areapower') & ~isempty(settings.areapower) & ...
        isfield(R,'areapower') & R.areapower > settings.areapower
    tmp = tmp + 1;
end
if isfield(settings,'bplast') & ~isempty(settings.bplast) & ...
        isfield(R,'bplast') & R.bplast > settings.bplast
    tmp = tmp + 1;
end
if isfield(settings,'percentbad') & ~isempty(settings.percentbad)
    if isfield(R,'badchans') & ~isempty(R.badchans) & length(R.badchans) > 1
        if isfield(settings,'REF') & isfield(settings.REF,'lap_percent') & ~isempty(settings.REF.lap_percent)
            R.bchans = sum(R.badchans > (settings.REF.lap_percent / 100));
        else
            R.bchans = sum(R.badchans>0);
        end
        if isempty(R.bchans) | R.bchans <= (length(R.badchans) * (settings.percentbad / 100))
            tmp = tmp+1;
        end
    else
        disp('    no info on bad channels, skip criteria')
        tmp = tmp + 1;
    end
else
    R.bchans = [];
end
if isfield(settings,'percentbadact') & ~isempty(settings.percentbadact)
    if isfield(R,'badact') & ~isempty(R.badact) & length(R.badact) > 1
        if R.badact(1,1) <= (R.badact(1,2)*(settings.percentbadact/100));
            tmp = tmp + 1;
        end
    else
        disp('    no info on bad activations, skip criteria')
        tmp = tmp + 1;
    end
end

if isfield(settings,'epochquality') & ~isempty(settings.epochquality)
    if isfield(R,'epochquality') & ~isempty(R.epochquality) & R.epochquality ~= 0
        if R.epochquality >= (settings.epochquality / 100)
            tmp = tmp + 1;
        end
    else
        disp('    no info on epoch quality, skip criteria')
        tmp = tmp + 1;
    end
end

if isfield(settings,'Npositiv') & ~isempty(settings.Npositiv) & tmp >= settings.Npositiv
    R.good = true;
else
    R.good = false;
end

end