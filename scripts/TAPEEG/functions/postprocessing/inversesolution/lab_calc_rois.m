function RIS = lab_calc_rois(data,header,ISmatrix,settings,ROIS,LOCS)
if ~isfield(settings,'roismethod')
    settings.roismethod = 'center';
end

% Prepare result matrix 
RIS.data = zeros(size(ROIS.solutionpts,2),size(data,2));
if isfield(ISmatrix,'x')
    RIS.x = zeros(size(ROIS.solutionpts,2),size(data,2));
    RIS.y = zeros(size(ROIS.solutionpts,2),size(data,2));
    RIS.z = zeros(size(ROIS.solutionpts,2),size(data,2));
end

% Calculation for method pseudo electrodes, G Bogaarts edit 2016
if strcmp(settings.roismethod,'pseudoelectrodes')
    ind=cell2mat(ROIS.solutionpts);
    ISmatrixPseudoE = lab_reduce_is(ISmatrix,ind);
    RIS = lab_calculate_is(data,header,ISmatrixPseudoE,settings.resultscalar);
    RIS.locs.x = LOCS.x(ind);
    RIS.locs.y = LOCS.y(ind);
    RIS.locs.z = LOCS.z(ind);
    RIS.locs.label = ROIS.label;
    RIS.locs.labels = ROIS.labels;
    RIS.locs = lab_locs2sph(RIS.locs);
    clearvars ISmatrixPseudoE ind  
end

% Calculation for method maxpower
if strcmp(settings.roismethod,'maxpower')
    for i = 1:size(ROIS.solutionpts,2)
        ISmatrixMaxP = lab_reduce_is(ISmatrix,ROIS.solutionpts{1,i});
        SPECT = lab_calculate_is(data,header,ISmatrixMaxP,2,1);
        if ~exist('n','var')
            n = find(SPECT.spectrum >= settings.spectralband(1,1), 1 );
            m = find(SPECT.spectrum < settings.spectralband(1,2), 1, 'last' );
        end
        [~,tmp] = max(mean(SPECT.data(:,n:m),2));
        maxsp(i) = ROIS.solutionpts{1,i}(1,tmp(1)); %#ok<AGROW>
        clearvars tmp
        slocs.x(1,i) = LOCS.x(1,maxsp(i));
        slocs.y(1,i) = LOCS.y(1,maxsp(i));
        slocs.z(1,i) = LOCS.z(1,maxsp(i));
        slocs.label{1,i} = ROIS.label{1,i};
        slocs.labels{1,i} = ROIS.labels{1,i};
        clearvars SPECT ISmatrixMaxP
    end
    ISmatrixMaxP = lab_reduce_is(ISmatrix,maxsp);
    RIS = lab_calculate_is(data,header,ISmatrixMaxP,settings.resultscalar);
    RIS.locs = slocs;
    clearvars ISmatrixMax maxsp m n slocs i
end

% Calculation for method mean
if strcmp(settings.roismethod,'mean')
    for i = 1:size(ROIS.solutionpts,2)
        % Reduce weigth matrix
        ISmatrixMean = lab_reduce_is(ISmatrix,ROIS.solutionpts{1,i});
        if isfield(ISmatrix,'x')
            if isfield(ISmatrixMean,'matrix')
                ISmatrixMean = rmfield(ISmatrix,'matrix');
            end
            ISmatrixMean.x = mean(ISmatrixMean.x,1);
            ISmatrixMean.y = mean(ISmatrixMean.y,1);
            ISmatrixMean.z = mean(ISmatrixMean.z,1);
        else
            disp('    warning: no control for optimal orientation possible')
            ISmatrixMean.matrix = mean(ISmatrixMean.matrix,1);
        end
        ISmatrixMean.TSolutionPointName = ROIS.label(1,i);
        ISmatrixMean.numsolutionpoints = 1;
        
        % Calculate IS
        RIStmp = lab_calculate_is(data,header,ISmatrixMean,settings.resultscalar);
        clearvars ISmatrixMean
        RIS.data(i,:) = RIStmp.data;
        RIS.locs.x(1,i) = mean(LOCS.x(1,ROIS.solutionpts{1,i}));
        RIS.locs.y(1,i) = mean(LOCS.y(1,ROIS.solutionpts{1,i}));
        RIS.locs.z(1,i) = mean(LOCS.z(1,ROIS.solutionpts{1,i}));
        RIS.locs.label{1,i} = ROIS.label{1,i};
        RIS.locs.labels{1,i} = ROIS.labels{1,i};
        if isfield(RIStmp,'x')
            RIS.x(i,:) = mean(RIStmp.x,1);
            RIS.y(i,:) = mean(RIStmp.y,1);
            RIS.z(i,:) = mean(RIStmp.z,1);
        end
    end
    clearvars i
end

% Calculation for method center
if strcmp(settings.roismethod,'center')
    for i = 1:size(ROIS.solutionpts,2)
        tmpx = LOCS.x(1,ROIS.solutionpts{1,i});
        tmpy = LOCS.y(1,ROIS.solutionpts{1,i});
        tmpz = LOCS.z(1,ROIS.solutionpts{1,i});
        for j = 1:size(tmpx,2)
            distance(1,j) = ((tmpx(1,j) - mean(tmpx))^2 + (tmpy(1,j) - mean(tmpy))^2 + (tmpz(1,j) - mean(tmpz))^2)^0.5; %#ok<AGROW>
        end
        tmp = find(distance == min(distance));
        sploc(i) = ROIS.solutionpts{1,i}(1,tmp(1,1)); %#ok<AGROW>
        clearvars tmpx tmpy tmpz distance tmp
        slocs.x(1,i) = LOCS.x(1,sploc(i));
        slocs.y(1,i) = LOCS.y(1,sploc(i));
        slocs.z(1,i) = LOCS.z(1,sploc(i));
        slocs.label{1,i} = ROIS.label{1,i};
        slocs.labels{1,i} = ROIS.labels{1,i};
    end
    ISmatrixCent = lab_reduce_is(ISmatrix,sploc);
    RIS = lab_calculate_is(data,header,ISmatrixCent,settings.resultscalar);
    RIS.locs = slocs;
    clearvars ISmatrixCent sploc slocs i
end
RIS.numsolutionpoints=size(RIS.data,1);
RIS.numtimeframes=size(RIS.data,2);
RIS.numchannels = size(RIS.data,1);
RIS.numdatachannels = size(RIS.data,1);
RIS.channels = char(RIS.locs.label);
RIS.labels = RIS.locs.labels;
RIS.numauxchannels = 0;