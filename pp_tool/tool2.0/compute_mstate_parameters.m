function [ output ] = compute_mstate_parameters( confstruct, eeg, maps)
tic
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%     eeg = avgref(eeg, confstruct.nch);
    dur_state = cell(size(maps,1), size(1: confstruct.tf : size(eeg,1),2));
    gfp_state = cell(size(maps,1), size(1: confstruct.tf : size(eeg,1),2));
%     freq_state = cell(size(maps,1), size(1: confstruct.tf : size(eeg,1),2));
    freq_state = [];
    mean_dur = [];
    time_cov = [];
    mean_gfp = [];
    for i = 1: confstruct.tf : size(eeg,1)
        if i+confstruct.tf-1 > length(eeg)
            epoch = eeg(i:end,:);
            ep_idx = round(i/confstruct.tf)+1;
        else
            ep_idx = round(i/confstruct.tf)+1;
            epoch = eeg(i:(i+confstruct.tf-1), :);
        end
        
        gfp_curve = compute_gfp(epoch, confstruct.method_GFPeak);
        if (confstruct.debug)
            disp(strcat('epoch ', num2str(ep_idx), ' gfp computed'))
        end       
        % detrend
        epoch_mean = mean(epoch,2);
        epoch = epoch - repmat(epoch_mean, 1, confstruct.nch);
        if (confstruct.debug)
            disp(strcat('epoch ', num2str(ep_idx), ' column mean subtracted'))
        end
        % detrend
        maps_mean = mean(maps, 2);
        maps = maps - repmat(maps_mean, 1, confstruct.nch);
        if (confstruct.debug)
            disp(strcat('epoch ', num2str(ep_idx), ' maps column mean subtracted'))
        end
        % set  gfp to 1 epoch
        epoch = set_gfp_all_1(epoch, gfp_curve);
        % set  gfp to 1 maps
        gfp_curve_maps = compute_gfp(maps, confstruct.method_GFPeak);
        maps = set_gfp_all_1(maps, gfp_curve_maps);
        if (confstruct.debug)
            disp(strcat('epoch ', num2str(ep_idx), ' gfp set to 1'))
        end
        % gfp peak
        [gfp_peaks_indices, gfp_peaks_values, gfp_curve] = compute_gfp_peaks(gfp_curve, confstruct.use_gfp_peaks, confstruct.maxima_method, confstruct.minDistance, confstruct.minAmplitude);
        if (confstruct.showpeaks)
            plot([1:length(gfp_curve)],gfp_curve,gfp_peaks_indices,gfp_peaks_values,'*')
        end
        if (confstruct.debug)
            disp(strcat('epoch ', num2str(ep_idx), ' gfp indices extracted'))
        end
        % compear each gfp eeg corresponding peak to maps
        map = containers.Map('KeyType','double','ValueType','any');
        map_index = containers.Map('KeyType','double','ValueType','any');
        map_corr = containers.Map('KeyType','double','ValueType','any');
        for p = 1:size(gfp_peaks_indices, 2)
            peak = gfp_peaks_indices(p);
            peak_tf = epoch(peak,:);
            for m = 1: size(maps,1)
                if (strcmp(confstruct.similarity_measure, 'correlation'))
                    r = corrcoef(maps(m,:),peak_tf);
%                     idx = setdiff([1,2,3,4], m);
%                     pr  = partialcorr([peak_tf',maps(m,:)'],maps(idx,:)');
                    single_gfpm_corr(m) = abs(r(1,2));
%                     single_gfpm_pcorr(m) = abs(pr(1,2));
                elseif (strcmp(confstruct.similarity_measure, 'dissimilarity'))
                    single_gfpm_diss(m) = globalDissimilarity(maps(m,:),peak_tf);
                end
            end
            map(peak) =  single_gfpm_corr;
            [val, ind] = max(single_gfpm_corr,[],2);
            map_index(peak) = ind;
            map_corr(peak) = val;
            single_epoch_corr(p,:)=single_gfpm_corr;
%             single_epoch_pcorr(p,:)=single_gfpm_pcorr;
        end
        [~,IND] = max(single_epoch_corr,[],2);
%         [~,INDp] = max(single_epoch_pcorr,[],2);
        previous = gfp_peaks_indices(1);
        start = gfp_peaks_indices(1);
        start_state = [start];
        corr_state = [map_corr(start)];
        ms_state = [map_index(start)];
        for p = 1: size(gfp_peaks_indices,2)
            p_tf = gfp_peaks_indices(p);
            if map_index(p_tf) == map_index(previous)
                previous = p_tf;
                continue 
            else
                start = previous + ((p_tf-previous)/2);
                start_state = [start_state; start];
                ms_state = [ms_state; map_index(p_tf)];
                corr_state = [corr_state; map_corr(p_tf)];
                previous = p_tf; 
            end      
        end
        %durata per ogni epoca
        k = 1;
        states_sequence = zeros(size(gfp_curve));
        states_sequence(1:start_state(1))=ms_state(1);
        states_sequence(start_state(end):length(gfp_curve))=ms_state(end);
        while k <= length(start_state)-1
            beg  = start_state(k);
            stop = start_state(k+1);
            dur = (stop - beg);
            state = ms_state(k);
            states_sequence(beg:stop)=state;
            if dur > confstruct.segment_min_len             
                mean_gfp = mean(gfp_curve(beg:stop));
                dur_state{state, ep_idx} = [dur_state{state, ep_idx}, dur*(1000/confstruct.fs)];
%                 dur_state{state, ep_idx} = [dur_state{state, ep_idx}, dur];
                gfp_state{state, ep_idx} = [gfp_state{state, ep_idx}, mean_gfp]; 
            else
%                 disp(strcat('epoca: ', num2str(ep_idx)))
%                 disp(strcat('durata: ',num2str(dur)))
%                 disp('REMOVED SMALL SEGMENT')
            end
            k=k+1;
        end
        %
        %durata totale
        dur_sum = 0;
        for ii = 1: size(dur_state,1)
            if ~isempty(dur_state{ii, ep_idx}) 
                dur_sum = dur_sum + sum(dur_state{ii, ep_idx});
            end
        end
        %frequenza
        for mapnr = 1: size(maps,1)
            if isempty(dur_state{mapnr, ep_idx})
                if (confstruct.debug)
                    disp(strcat('In epoch ', num2str(ep_idx), ' maps ', num2str(mapnr), ' does not occur. Frequency set to 0.' ));
                end
                freq_state(mapnr, ep_idx) = 0;
            else
%                 freq_state(mapnr, ep_idx) = length(dur_state{mapnr, ep_idx})*(confstruct.tf*(1000/confstruct.fs/dur_sum)/(confstruct.tf/confstruct.fs));
%                 freq_state(mapnr, ep_idx) = (length(dur_state{mapnr, ep_idx})*confstruct.fs)/confstruct.tf;
                freq_state(mapnr, ep_idx) = (length(dur_state{mapnr, ep_idx}));

            end
        end
        %mean durata, time coverage and gfp medio per ogni epoca
        for k = 1:size(dur_state,1)
            for j = 1:size(dur_state, 2)
                mean_dur(k, j) = nanmean(dur_state{k, j});
                time_cov(k, j) = sum(dur_state{k, j})/dur_sum;
                mean_gfp(k, j) = mean(gfp_state{k, j});
            end
        end
        nr_of_gfp(ep_idx) = length(gfp_peaks_indices);
        total_mean_of_gfp(ep_idx) = mean(gfp_curve);
        % varianza spiegata
        b_assignment = cell2mat(values(map_index));
        b_loading = cell2mat(values(map_corr));

        for m = 1: size(maps, 1)
            exp_var(m) = sum(b_loading(b_assignment==m)) / sum(std(epoch(gfp_peaks_indices,:),1, 2));
        end
        exp_var_all = sum(exp_var);
        GEV(ep_idx,:)=exp_var;
        GEV_ALL(1,ep_idx)=exp_var_all;
        if (confstruct.debug)
                disp(strcat('epoch ', num2str(ep_idx), ' all parameteres computed'))
        end
        
        clear dur 
        
    end
   
    output.mean_dur = mean_dur;
    output.tot_dur = dur_sum;
    output.freq = freq_state;
    output.gev = exp_var;
    output.GEV_ALL = GEV_ALL;
    output.cov = time_cov;
    output.nr_of_gfp_points = nr_of_gfp;
    output.mean_gfp = total_mean_of_gfp;
    output.states_sequence=states_sequence;
    output.single_epoch_corr=single_epoch_corr;
%     output.single_epoch_pcorr=single_epoch_pcorr;
    
%     clear mean_dur start_state 
toc 
end






