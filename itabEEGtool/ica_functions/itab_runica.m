function [ EEG ] = itab_runica( EEG )
%ITAB_ICA Summary of this function goes here
%   Detailed explanation goes here
    sigs = EEG.ICA.icadata;
    [Weigth, Sphere] = runica(sigs);
    W = Weigth*Sphere;
    A = inv(W);
    IC = W*sigs;
    [chan, Nc] = size(A);
    [Nc,len] = size(IC);
    if Nc > 0
        for i=1:Nc
            potenza(i) = mean(A(:,i).^2);
        end
        [potenza, order] = sort(potenza);
        potenza = fliplr(potenza);
        order = fliplr(order);
        A = A(:, order);
        W = W(order, :);
        IC = IC(order, :);
    end
    EEG.ICA.A = A;
    EEG.ICA.W = W;
    EEG.ICA.IC = IC;
end

