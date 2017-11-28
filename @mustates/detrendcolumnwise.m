function [eeg] = detrendcolumnwise(self, eeg)
for i  = 1:self.confstruct.tf:size(eeg,1)
     if i+ self.confstruct.tf-1 > size(eeg)
         epoch = eeg(i:end,:);
         epoch_mean_subtracted = epoch - repmat(mean(epoch), size(epoch,1), 1);
         eeg(i:end, :) = epoch_mean_subtracted;
     else
         epoch = eeg(i:(i + self.confstruct.tf-1), :);
         epoch_mean_subtracted = epoch - repmat(mean(epoch), self.confstruct.tf, 1);
         eeg(i:(i + self.confstruct.tf-1), :) = epoch_mean_subtracted;
     end        
end