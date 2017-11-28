function [MSClass,MSFit,gfp] = FitMicrostatesCorr(data,MSMaps,Smooth,b,lambda,Start,Stop)

        % output arguments
        % MSClass = assigned time frame to a max correlate template;
        % MSFi = value of correlation of the assigned time frame;
        % gfp = gfp of each assgned time frame

        % input arguments
        % data = data to be fitted (n_channels*n_time_frames);
        % MSMaps = templates (maps) to fit;
        % Smooth = temporal smoothing: 1 = yes; 0 = no;
        % Start = initial time frame of the segment of data to be fitted;
        % Stop = final time frame of the segment of data to be fitted;

        MSMaps = NormDimL2(MSMaps,2) / sqrt(size(MSMaps,2));

        ne = size(data,1);
        nt = size(data,2);

        gfp = std(data,1);

        MSFit   = zeros(1,nt);
        MSClass = zeros(1,nt);

        if Smooth == 1 && size(MSMaps,1) > 1
            [idx,bf] = RaguSmoothLabels(data,MSMaps,b,lambda,Start,Stop );
        else
%             [bfall,idxall] = max(MSMaps * data,[],1);
            for m = 1:size(MSMaps,1)
                for i = 1:size(data, 2)
                     CC = corrcoef(MSMaps(m,:),data(:,i));
                     CORR(m,i)=abs(CC(1,2));
                end
            end
            [bfall,idxall] = max(CORR,[],1);
            bf = zeros(1,size(data,2));
            idx = bf;
            bf(Start:Stop) = bfall(Start:Stop);
            idx(Start:Stop) = idxall(Start:Stop);
        end
        MSClass(1,:) = idx;
        MSFit(1,:)   = bf;

MSFit = MSFit / sqrt(ne);
