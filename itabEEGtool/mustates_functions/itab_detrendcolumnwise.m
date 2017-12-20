function [eeg] = detrendcolumnwise(eeg, confstruct)
for i  = 1:confstruct.tf:size(eeg,1)
     if i+confstruct.tf-1 > size(eeg)
         epoch = eeg(i:end,:);
         epoch_mean_subtracted = epoch - repmat(mean(epoch), size(epoch,1), 1);
         eeg(i:end, :) = epoch_mean_subtracted;
     else
         epoch = eeg(i:(i + confstruct.tf-1), :);
         epoch_mean_subtracted = epoch - repmat(mean(epoch), confstruct.tf, 1);
         eeg(i:(i + confstruct.tf-1), :) = epoch_mean_subtracted;
     end        
end