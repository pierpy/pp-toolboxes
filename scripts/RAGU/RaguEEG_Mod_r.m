function[b_model,b_ind,b_loading,exp_var] = RaguEEG_Mod_r(eeg,n_mod,reruns,max_n,flags)

%function[b_model,b_ind,b_loading,exp_var] = RaguEEG_Mod_r(eeg,n_mod,reruns,max_n,flags)

% EEG_MOD Create the EEG model I'm working with
%
% function[b_model,b_ind,b_loading,exp_var] = eeg_mod_r(eeg,n_mod,reruns,max_n)

% input arguments
% eeg = the input data (number of time instances * number of channels)
% n_mod = the number of microstate clusters that you want to extract
% reruns = the number of reiterations (use about 20)
% max_n = maximum number of eeg timepoints used for cluster indentification

% output arguments
% b_model = cluster centers (microstate topographies)
% b_ind = cluster assignment for each moment of time
% b_loading = Amplitude of the assigned cluster at each moment in time
% exp_var = explained variance of the model

% Copyright 2009-2011 Thomas Koenig
% distributed under the terms of the GNU AFFERO General Public License

if (size(n_mod,1) ~= 1)
	error('Second argument must be a scalar')
end

if (size(n_mod,2) ~= 1)
	error('Second argument must be a scalar')
end

[n_frame,n_chan] = size(eeg);
h = eye(n_chan)-1/n_chan;
eeg = eeg*h;									% Average reference of data 

if nargin < 3
    reruns = 1;
end

if nargin < 4
    max_n = n_frame;
end

if isempty(max_n)
    max_n = n_frame;
end

if isempty(findstr(flags,'p'))
    pmode = 0;
else
    pmode = 1;
end

org_data = eeg;
best_fit = 0;

tic

hndl = waitbar(0,sprintf('Fitting %i microstates, please wait...',n_mod));


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
        if pmode
            covm    = eeg * model';						% Get the unsigned covariance matrix
        else
            covm    = abs(eeg * model');						% Get the unsigned covariance matrix
        end
        [c,ind] =  max(covm,[],2);				            % Look for the best fit

        for i = 1:n_mod
            idx = find (ind == i);
            if pmode
                model(i,:) = mean(eeg(idx,:));
            else
                cvm = eeg(idx,:)' * eeg(idx,:);
                [v] = pc_evectors(cvm,1);
                model(i,:) = v(:,1)';
            end
        end
		model   = normr(model)*h;						% Average Reference, equal variance of model
        covm    = eeg*model';							% Get the unsigned covariance 
        if pmode
            [c,ind] =  max(covm,[],2);				% Look for the best fit
        else
            [c,ind] =  max(abs(covm),[],2);				% Look for the best fit
        end
    end % while any
    covm    = org_data*model';							% Get the unsigned covariance 
    if pmode
        [loading,ind] =  max(covm,[],2);				% Look for the best fit
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
    
    waitbar(run/reruns,hndl);
    set(hndl,'Name',sprintf('Remaining time: %01.0f:%02.0f min',floor(toc()*(reruns/run)/60),rem(toc()*(reruns/run-1),60)));
    
end % for run = 1:reruns

close(hndl);
