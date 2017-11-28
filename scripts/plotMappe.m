for k=1:size(xx)
    subplot(2,2,k)
    topoplot(xx(k,:), EEG.chanlocs)
end