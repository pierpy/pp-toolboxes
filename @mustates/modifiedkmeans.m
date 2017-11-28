function [b_model, b_ind, b_loading, exp_var] = modifiedkmeans(self, eeg, n_mod, reruns)
    % output arguments
        % b_model = cluster centers (microstate topographies)
        % b_ind = cluster assignment for each moment of time
        % b_loading = Amplitude of the assigned cluster at each moment in time
        % exp_var = explained variance of the model

    % input arguments
        % eeg = data to be analyzed (n_time_frames*n_channels)
        % n_mod = number of cluster with wich decompose the data
        % self.confstruct.pmode = polarity mode: 0 polarity indipendent; 1 polarity dependent
        % average = average reference: 0 no; 1 yes.

    if (size(n_mod,1) ~= 1)
        error('Second argument must be a scalar')
    end

    if (size(n_mod,2) ~= 1)
        error('Second argument must be a scalar')
    end

    [n_frame,n_chan] = size(eeg);
    h = eye(n_chan)-1/n_chan;
    eeg = eeg*h;	

    max_n=n_frame;
    org_data = eeg;
    best_fit = 0;

    for run = 1:reruns
        if nargin > 3
            idx = randperm(n_frame);
            eeg = org_data(idx(1:max_n),:);
        end

        idx = randperm(max_n);
        model = eeg(idx(1:n_mod),:);
        model   = normr(model)*h;							% Average Reference, equal variance of model

        o_ind   = zeros(max_n,1);							% Some initialization
        ind     =  ones(max_n,1);
        count   = 0;

        while any(o_ind - ind)
            count   = count+1;
            o_ind   = ind;
            if self.confstruct.pmode
                covm    = eeg * model';						% Get the unsigned covariance matrix
            else
                covm    = abs(eeg * model');						% Get the unsigned covariance matrix
            end
            [c,ind] =  max(covm,[],2);				            % Look for the best fit

            for i = 1:n_mod
                idx = find (ind == i);
                if self.confstruct.pmode
                    model(i,:) = mean(eeg(idx,:));
                else
                    cvm = eeg(idx,:)' * eeg(idx,:);
                    [v] = pc_evectors(cvm,1,0);
                    model(i,:) = v(:,1)';
                end
            end

            model   = normr(model)*h;	% Average Reference, equal variance of model
            covm    = eeg*model';							% Get the unsigned covariance 
            if self.confstruct.pmode
                  [c,ind] =  max(covm,[],2);	
     %            [c,ind] =  max(abs(covm),[],2);	% Look for the best fit
            else
                  [c,ind] =  max(abs(covm),[],2);				% Look for the best fit
            end
        end % while any
        covm    = org_data*model';							% Get the unsigned covariance 
        if self.confstruct.pmode
                [loading,ind] =  max(covm,[],2);	
    %           [loading,ind] =  max(abs(covm),[],2);	% Look for the best fit
        else
                [loading,ind] =  max(abs(covm),[],2);				% Look for the best fit
        end

        tot_fit = sum(loading);
        if (tot_fit > best_fit)
            b_model   = model;
            b_ind     = ind;
            b_loading = loading/sqrt(n_chan);
            best_fit  = tot_fit;
            exp_var = sum(b_loading)/sum(std(eeg,1,2));
        end

    end
end% for run = 1:reruns


