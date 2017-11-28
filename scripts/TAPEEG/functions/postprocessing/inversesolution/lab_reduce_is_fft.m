% Reduce IS_FFT_Result (lab_inversesolution_fft) to ROIS
%
% written by F. Hatz 2012

function ISresult = lab_reduce_is_fft(ISresult,ROIS,cfg)
if ~isfield(cfg,'roismethod')
    cfg.roismethod = 'mean';
end

% Prepare result matrix
tmp.data = zeros(size(ROIS.solutionpts,2),size(ISresult.data,2));
if isfield(ISresult,'x')
    tmp.x = zeros(size(ROIS.solutionpts,2),size(ISresult.data,2));
    tmp.y = zeros(size(ROIS.solutionpts,2),size(ISresult.data,2));
    tmp.z = zeros(size(ROIS.solutionpts,2),size(ISresult.data,2));
end

% Calculation for method maxpower
if strcmp(cfg.roismethod,'maxpower')
    if isfield(ISresult,'spectrum')
        specdata = ISresult.data;
        for i = 1:size(ROIS.solutionpts,2)
            for j = 1:size(specdata,2)
                [~,maxsp] = max(specdata(ROIS.solutionpts{1,i},j));
                maxsp = ROIS.solutionpts{1,i}(1,maxsp);
                maxsp = maxsp(1,1);
                tmp.data(i,j) = ISresult.data(maxsp,j);
                if isfield(ISresult,'x')
                    tmp.x(i,j) = ISresult.x(maxsp,j);
                    tmp.y(i,j) = ISresult.y(maxsp,j);
                    tmp.z(i,j) = ISresult.z(maxsp,j);
                end
            end
            tmp.locs.x(1,i) = ISresult.locs.x(1,maxsp);
            tmp.locs.y(1,i) = ISresult.locs.y(1,maxsp);
            tmp.locs.z(1,i) = ISresult.locs.z(1,maxsp);
            tmp.locs.label{1,i} = ROIS.label{1,i};
            tmp.locs.labels{1,i} = ROIS.labels{1,i};
        end
    else
        disp('   MaxPower-Calculation not possible, using mean-method')
        cfg.roismethod = 'mean';
    end
end

% Calculation for method mean
if strcmp(cfg.roismethod,'mean')
    for i = 1:size(ROIS.solutionpts,2)
        tmp.data(i,:) = mean(ISresult.data(ROIS.solutionpts{1,i},:).^0.5).^2;
        if isfield(ISresult,'x')
            tmp.x(i,:) = mean(ISresult.x(ROIS.solutionpts{1,i},:).^0.5).^2;
            tmp.y(i,:) = mean(ISresult.y(ROIS.solutionpts{1,i},:).^0.5).^2;
            tmp.z(i,:) = mean(ISresult.z(ROIS.solutionpts{1,i},:).^0.5).^2;
        end
        tmp.locs.x(1,i) = mean(ISresult.locs.x(1,ROIS.solutionpts{1,i}));
        tmp.locs.y(1,i) = mean(ISresult.locs.y(1,ROIS.solutionpts{1,i}));
        tmp.locs.z(1,i) = mean(ISresult.locs.z(1,ROIS.solutionpts{1,i}));
        tmp.locs.label{1,i} = ROIS.label{1,i};
        tmp.locs.labels{1,i} = ROIS.labels{1,i};
    end
end

% Calculation for method center
if strcmp(cfg.roismethod,'center')
    for i = 1:size(ROIS.solutionpts,2)
        tmpx = ISresult.locs.x(1,ROIS.solutionpts{1,i});
        tmpy = ISresult.locs.y(1,ROIS.solutionpts{1,i});
        tmpz = ISresult.locs.z(1,ROIS.solutionpts{1,i});
        for j = 1:size(tmpx,2)
            distance(1,j) = ((tmpx(1,j) - mean(tmpx))^2 + (tmpy(1,j) - mean(tmpy))^2 + (tmpz(1,j) - mean(tmpz))^2)^0.5;
        end
        sploc = (distance == min(distance));
        sploc = ROIS.solutionpts{1,i}(1,sploc);
        clearvars tmpx tmpy tmpz distance
        sploc = sploc(1,1);
        tmp.data(i,:) = ISresult.data(sploc,:);
        tmp.locs.x(1,i) = ISresult.locs.x(1,sploc);
        tmp.locs.y(1,i) = ISresult.locs.y(1,sploc);
        tmp.locs.z(1,i) = ISresult.locs.z(1,sploc);
        tmp.locs.label{1,i} = ROIS.label{1,i};
        tmp.locs.labels{1,i} = ROIS.labels{1,i};
        if isfield(ISresult,'x')
            tmp.x(i,:) = ISresult.x(sploc,:);
            tmp.y(i,:) = ISresult.y(sploc,:);
            tmp.z(i,:) = ISresult.z(sploc,:);
        end
    end
end

ISresult.numsolutionpoints=size(tmp.data,1);
ISresult.numtimeframes=size(tmp.data,2);
ISresult.data = tmp.data;
if isfield(tmp,'x')
    ISresult.x = tmp.x;
    ISresult.y = tmp.y;
    ISresult.z = tmp.z;
end
ISresult.locs = tmp.locs;
ISresult.numchannels = size(tmp.data,1);
ISresult.numdatachannels = size(tmp.data,1);
ISresult.channels = char(tmp.locs.label);
ISresult.labels = tmp.locs.labels;
ISresult.numauxchannels = 0;
clearvars tmp