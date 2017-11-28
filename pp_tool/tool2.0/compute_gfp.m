function gfp_curve = compute_gfp(eeg, method)
    ntf = size(eeg,1);
    gfp_curve = zeros(1,ntf);
    if strcmp(method,'GFPL2')
        for i = 1 : ntf
            x = eeg(i,:);
            gfp = sqrt(sum((x-mean(x)).^2)/length(x));
            gfp_curve(i)=gfp;
        end
    elseif strcmp(method,'GFP1')
        for i = 1 : ntf
            x = eeg(i,:);
            gfp = sum(abs(x-mean(x)))/length(x);
            gfp_curve(i)=gfp;
        end
    elseif strcmp(method, 'STD')
        for i = 1 : ntf
            x = eeg(i,:);
            gfp = std(x);
            gfp_curve(i)=gfp;
        end
    end