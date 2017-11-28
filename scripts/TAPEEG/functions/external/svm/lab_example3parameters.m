function [C, sigma] = lab_example3parameters(data,result,cfg)
% Initialization
C = 1;
sigma = 0.3;
% This function returns the optimal C and sigma learning parameters 
% found using the cross validation set. I have used svmPredict to predict 
% the labels on the cross validation set. For example, predictions = svmPredict(model, Xval);
% will return the predictions on the cross validation set.
% Prediction error can be computed by using mean(double(predictions ~= yval))

ntrials = size(data,1);
error = zeros(8,8);
set = [0.01,0.03,0.1,0.3,1,3,10,30];
for nloop = 1:cfg.SVM.loops
    randvector = randperm(ntrials);
    cut = ceil(ntrials*(1-cfg.SVM.percenttest));
    X = data(randvector(1:cut),:);
    y = result(randvector(1:cut),1);
    if cut < size(data,1)
        Xval = data(randvector(cut+1:end),:);
        yval = result(randvector(cut+1:end),1);
    else
        return
    end
    for i= 1:8
        c(i) = set(i);
        for j = 1:8
            ss(j) = set(j);
            
            model = svmTrain(X, y, c(i), @(x1, x2) gaussianKernel(x1, x2, ss(j)));
            
            predictions  = svmPredict(model,Xval);
            
            error(i,j) = error(i,j) + mean(double(predictions~=yval));
        end
    end
end

%error;
x = min(min(error));
for i =1:8
    for j =1:8
        if error(i,j) == x
            C = set(i);
            sigma = set(j);
            break;
        end
    end
end

end
