function eeg = substract_column_mean_epochwise(eeg, tf)
    for i  = 1:tf:size(eeg,1)
         if i+tf-1 > size(eeg)
             epoch = eeg(i:end,:);
             epoch_mean_subtracted = epoch - repmat(mean(epoch), size(epoch,1), 1);
             eeg(i:end, :)=epoch_mean_subtracted;
         else
             epoch = eeg(i:(i+tf-1), :);
             epoch_mean_subtracted = epoch - repmat(mean(epoch), tf, 1);
             eeg(i:(i+tf-1), :)=epoch_mean_subtracted;
         end        
    end